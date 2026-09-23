using Serialization, Statistics
using JLD2

include("Validation.jl")
using .Validation: write_rows

const CONFIG_PATTERN =
    r"^cfg_(\d+)_(fixed|inverse_sqrt|cos)_i(\d+)_lr(0p5|1p0|2p0)$"

finite_value(value) =
    !ismissing(value) && value isa Number && isfinite(value)

function config_metadata(path, root)
    parts = splitpath(relpath(dirname(path), root))
    config_position = findfirst(part -> match(CONFIG_PATTERN, part) !== nothing,
                                parts)
    isnothing(config_position) &&
        error("cannot find parameter configuration in path: $path")
    matched = match(CONFIG_PATTERN, parts[config_position])
    phase = Symbol(parts[1])
    phase in (:method1, :method2, :reduction) ||
        error("unknown parameter-study phase in $path")
    (
        phase=String(phase),
        config=parse(Int, matched.captures[1]),
        schedule=matched.captures[2],
        iterations=parse(Int, matched.captures[3]),
        learning_rate_scale=parse(Float64,
            replace(matched.captures[4], "p" => ".")),
    )
end

function numeric_summary(rows, property)
    values = Float64[]
    for row in rows
        hasproperty(row, property) || continue
        value = getproperty(row, property)
        finite_value(value) && push!(values, Float64(value))
    end
    isempty(values) ?
        (minimum=missing, median=missing, mean=missing, maximum=missing) :
        (minimum=minimum(values), median=median(values),
         mean=mean(values), maximum=maximum(values))
end

function summarize(grouped)
    output = NamedTuple[]
    for (key, rows) in grouped
        residual = numeric_summary(rows, :residual)
        objective = numeric_summary(rows, :objective)
        relaxed_residual = numeric_summary(rows, :relaxed_residual)
        seconds = numeric_summary(rows, :seconds)
        statuses = hasproperty(first(rows), :status) ?
                   getproperty.(rows, :status) : fill("", length(rows))
        elements = numeric_summary(rows, :elements)
        push!(output, (
            phase=key.phase, case_name=key.case_name,
            config=key.config, schedule=key.schedule,
            iterations=key.iterations,
            learning_rate_scale=key.learning_rate_scale,
            attempts=length(rows),
            converged=count(==("converged"), statuses),
            not_converged=count(==("not_converged"), statuses),
            inadmissible=count(==("inadmissible"), statuses),
            failed=count(==("failed"), statuses),
            best_objective=objective.minimum,
            median_objective=objective.median,
            mean_objective=objective.mean,
            best_residual=residual.minimum,
            median_residual=residual.median,
            mean_residual=residual.mean,
            median_relaxed_residual=relaxed_residual.median,
            median_elements=elements.median,
            total_seconds=isempty(rows) ? missing :
                sum(Float64(getproperty(row, :seconds)) for row in rows
                    if hasproperty(row, :seconds) &&
                       finite_value(getproperty(row, :seconds))),
            median_seconds=seconds.median,
        ))
    end
    sort!(output; by=row -> (row.phase, row.config, row.case_name))
end

function parse_list(value, convert_value)
    ismissing(value) && return Any[]
    text = String(value)
    isempty(text) && return Any[]
    map(part -> isempty(part) ? missing : convert_value(part),
        split(text, ';'; keepempty=true))
end

function pareto_rows(rows)
    output = NamedTuple[]
    for row in rows
        hasproperty(row, :reduction_pareto) || continue
        flags = parse_list(row.reduction_pareto,
                           value -> lowercase(value) == "true")
        isempty(flags) && continue
        topologies = parse_list(row.reduction_topologies, identity)
        elements = parse_list(row.reduction_elements, value -> parse(Int, value))
        objectives = parse_list(row.reduction_objectives,
                                value -> parse(Float64, value))
        residuals = parse_list(row.reduction_residuals,
                               value -> parse(Float64, value))
        stiffness = parse_list(row.reduction_stiffness_errors,
                               value -> parse(Float64, value))
        scores = parse_list(row.reduction_selection_scores,
                            value -> parse(Float64, value))
        n = minimum(length.((flags, topologies, elements, objectives,
                            residuals, stiffness, scores)))
        for index in 1:n
            flags[index] === true || continue
            push!(output, (
                phase=row.phase, case_name=row.case_name,
                config=row.config, schedule=row.schedule,
                iterations=row.iterations,
                learning_rate_scale=row.learning_rate_scale,
                seed=row.seed, initialization=row.initialization,
                reduction_step=index - 1,
                selected=!ismissing(row.selected_reduction_step) &&
                         row.selected_reduction_step == index - 1,
                topology=topologies[index], elements=elements[index],
                objective=objectives[index], residual=residuals[index],
                stiffness_error=stiffness[index],
                selection_score=scores[index],
            ))
        end
    end
    output
end

function scalar_jld_metadata(path)
    JLD2.jldopen(path, "r") do file
        value(name, default=missing) =
            haskey(file, name) ? read(file, name) : default
        (
            converged=Bool(value("converged", false)),
            objective=value("objective"),
            optimization_objective=value("optimization_objective"),
            residual=value("residual"),
            method=String(value("method", "")),
            case_name=String(value("case_name", "")),
        )
    end
end

rank(metadata) = (
    metadata.converged ? 0 : 1,
    finite_value(metadata.objective) ? Float64(metadata.objective) : Inf,
    finite_value(metadata.residual) ? Float64(metadata.residual) : Inf,
    finite_value(metadata.optimization_objective) ?
        Float64(metadata.optimization_objective) : Inf,
)

root = abspath(length(ARGS) >= 1 ? ARGS[1] :
    joinpath(@__DIR__, "results", "parameter_study"))
output = abspath(length(ARGS) >= 2 ? ARGS[2] : joinpath(root, "aggregated"))
isdir(root) || error("parameter-study root does not exist: $root")

row_files = String[]
solution_files = String[]
for (directory, _, files) in walkdir(root), name in files
    path = joinpath(directory, name)
    name == "rows.jls" && push!(row_files, path)
    endswith(name, "_best_solution.jld2") && push!(solution_files, path)
end
isempty(row_files) && error("no rows.jls files found below $root")

all_rows = Dict{String,Vector{NamedTuple}}(
    "method1" => NamedTuple[],
    "method2" => NamedTuple[],
    "reduction" => NamedTuple[],
)
grouped = Dict{NamedTuple,Vector{NamedTuple}}()
observed_shards = Dict{NamedTuple,Set{Int}}()
expected_shards = Dict{NamedTuple,Int}()

for file in sort(row_files)
    payload = deserialize(file)
    metadata = config_metadata(file, root)
    key = (; metadata..., case_name=String(payload.case_name),
           method=String(payload.method))
    append!(all_rows[metadata.phase],
        [(; metadata..., case_name=String(payload.case_name),
           shard_index=payload.shard_index, row...) for row in payload.rows])
    append!(get!(grouped, key, NamedTuple[]),
            [(; metadata..., case_name=String(payload.case_name), row...)
             for row in payload.rows])
    push!(get!(observed_shards, key, Set{Int}()), payload.shard_index)
    expected_shards[key] = payload.shard_count
end

completeness = NamedTuple[]
for key in sort(collect(keys(expected_shards));
                by=key -> (key.phase, key.config, key.case_name))
    expected = expected_shards[key]
    present = observed_shards[key]
    missing_shards = setdiff(1:expected, present)
    push!(completeness, (; key..., expected_shards=expected,
        present_shards=length(present),
        complete=isempty(missing_shards),
        missing_shards=join(missing_shards, ";")))
end

mkpath(output)
for phase in ("method1", "method2", "reduction")
    rows = all_rows[phase]
    isempty(rows) || write_rows(joinpath(output, "$(phase)_runs.csv"), rows)
end
summary = summarize(grouped)
write_rows(joinpath(output, "parameter_study_summary.csv"), summary)
write_rows(joinpath(output, "parameter_study_completeness.csv"), completeness)

reduction_pareto = pareto_rows(all_rows["reduction"])
isempty(reduction_pareto) ||
    write_rows(joinpath(output, "parameter_study_pareto.csv"),
               reduction_pareto)

best_directory = joinpath(output, "best_solutions")
mkpath(best_directory)
best = Dict{String,Tuple{Any,String}}()
for file in solution_files
    metadata = config_metadata(file, root)
    candidate = try
        scalar_jld_metadata(file)
    catch err
        @warn "Skipping unreadable solution" file exception=(err, catch_backtrace())
        continue
    end
    phase = metadata.phase
    if !haskey(best, phase) || rank(candidate) < rank(first(best[phase]))
        best[phase] = (candidate, file)
    end
end
for (phase, (metadata, source)) in best
    destination = joinpath(best_directory,
        "$(phase)_$(metadata.case_name)_best_solution.jld2")
    cp(source, destination; force=true)
end

incomplete = count(row -> !row.complete, completeness)
println("Aggregated $(length(row_files)) shard files into $output")
println("Configurations: $(length(completeness)); incomplete: $incomplete")
println("Raw shard files were not modified.")
incomplete == 0 || exit(2)
