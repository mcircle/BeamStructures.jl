using Serialization, TOML

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
for (directory, _, names) in walkdir(root), name in names
    name == "rows.jls" && push!(files, joinpath(directory, name))
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

println("Merged $(length(files)) shards into $output")
