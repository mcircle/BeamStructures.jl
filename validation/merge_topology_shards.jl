using Serialization, TOML, JLD2

include("Validation.jl")
include("topology_generation.jl")
include("topology_evaluation.jl")

using .Validation
using .TopologyGeneration
using .TopologyEvaluation

root = isempty(ARGS) ? joinpath(@__DIR__, "results", "lsf_shards") :
       abspath(ARGS[1])
output = length(ARGS) >= 2 ? abspath(ARGS[2]) :
         joinpath(@__DIR__, "results", "lsf_merged")
isdir(root) || error("shard directory not found: $root")
mkpath(output)

settings = TOML.parsefile(joinpath(@__DIR__, "config.toml"))
include("topology_adapter.jl")
cases = topology_study_cases(settings)
case_lookup = Dict(case.name => case for case in cases)

files = String[]
candidate_files = String[]
for (directory, _, names) in walkdir(root), name in names
    name == "rows.jls" && push!(files, joinpath(directory, name))
    endswith(name, "_best_solution.jld2") &&
        push!(candidate_files, joinpath(directory, name))
end
isempty(files) && error("no rows.jls shard files found below $root")

groups = Dict{Tuple{Symbol,String},Vector{Any}}()
for file in files
    payload = deserialize(file)
    payload.schema_version == 1 || error("unsupported shard schema in $file")
    key = (payload.method, payload.case_name)
    push!(get!(groups, key, Any[]), payload)
end

merged = Dict{Tuple{Symbol,String},Vector{NamedTuple}}()
for (key, shards) in groups
    counts = unique(getproperty.(shards, :shard_count))
    length(counts) == 1 || error("inconsistent shard counts for $key")
    expected = only(counts)
    indices = sort(getproperty.(shards, :shard_index))
    indices == collect(1:expected) ||
        error("missing or duplicate shards for $key: found $indices")
    rows = NamedTuple[]
    for shard in sort(shards; by=x -> x.shard_index)
        append!(rows, shard.rows)
    end
    if key[1] == :method1
        sort!(rows; by=row -> (row.topology, row.seed))
    else
        sort!(rows; by=row -> row.seed)
    end
    merged[key] = rows
end

for case in cases, method in (:method1, :method2)
    haskey(merged, (method, case.name)) ||
        error("$method shards missing for $(case.name)")
end

for name in readdir(output)
    path = joinpath(output, name)
    isfile(path) || continue
    if endswith(name, ".csv") || endswith(name, ".jld2") ||
       name in ("metadata.toml", "Manifest.toml")
        rm(path)
    end
end

topologies = enumerate_topologies()
points = settings["evaluation_points"]
write_study_inputs(topologies, cases, points, output)
record_environment(output; settings)

summaries = Dict{String,Vector{NamedTuple}}()
residual_limit = get(settings, "topology_residual_limit", Inf)
for case in cases
    name = case.name
    method1_key = (:method1, name)
    method2_key = (:method2, name)
    haskey(merged, method1_key) ||
        error("method1 shards missing for $name")
    haskey(merged, method2_key) ||
        error("method2 shards missing for $name")

    method1 = merged[method1_key]
    method2 = merged[method2_key]
    write_rows(joinpath(output, "$(name)_method1_runs.csv"), method1)
    write_rows(joinpath(output, "$(name)_method2_runs.csv"), method2)
    summary = summarize_topologies(method1; residual_limit)
    summaries[name] = summary
    write_rows(joinpath(output, "$(name)_topology_summary.csv"), summary)
    compare_method2(summary, method2; directory=output, name)
end

function stored_rank(data)
    converged = get(data, "converged", false)
    objective = get(data, "optimization_objective", Inf)
    residual = get(data, "residual", Inf)
    (converged ? 0 : 1,
     ismissing(objective) || !isfinite(objective) ? Inf : Float64(objective),
     ismissing(residual) || !isfinite(residual) ? Inf : Float64(residual))
end

best_files = Dict{Tuple{Symbol,String},Tuple{Any,String}}()
for file in candidate_files
    data = JLD2.load(file)
    key = (Symbol(data["method"]), String(data["case_name"]))
    haskey(groups, key) || continue
    if !haskey(best_files, key) ||
       stored_rank(data) < stored_rank(first(best_files[key]))
        best_files[key] = (data, file)
    end
end
for ((method, case_name), (_, source)) in best_files
    destination = joinpath(output,
        "$(case_name)_$(method)_best_solution.jld2")
    cp(source, destination; force=true)
end

println("Merged $(length(files)) shards into $output")

cleanup = lowercase(get(ENV, "BEAM_STUDY_CLEAN_SHARDS", "true")) in
          ("1", "true", "yes")
if cleanup
    root_path = realpath(root)
    output_path = abspath(output)
    output_below_root = startswith(output_path,
        root_path * string(Base.Filesystem.path_separator))
    protected = (root_path == output_path || root_path == homedir() ||
                 root_path == pwd() || dirname(root_path) == root_path ||
                 output_below_root)
    protected && error("refusing to remove unsafe shard directory: $root_path")
    rm(root_path; recursive=true)
    println("Removed merged shard directory $root_path")
end
