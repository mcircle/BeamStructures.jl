using Serialization, TOML

include("Validation.jl")
include("topology_generation.jl")
include("topology_evaluation.jl")

using .Validation
using .TopologyGeneration
using .TopologyEvaluation

required(name) = haskey(ENV, name) ? ENV[name] :
                 error("required environment variable $name is missing")

method = Symbol(required("BEAM_STUDY_METHOD"))
method in (:method1, :method2) ||
    error("BEAM_STUDY_METHOD must be method1 or method2")
case_name = required("BEAM_STUDY_CASE")
shard_index = parse(Int, required("BEAM_STUDY_SHARD_INDEX"))
shard_count = parse(Int, required("BEAM_STUDY_SHARD_COUNT"))
1 <= shard_index <= shard_count ||
    error("shard index must lie in 1:shard_count")

settings = TOML.parsefile(joinpath(@__DIR__, "config.toml"))
include("topology_adapter.jl")
cases = Dict(case.name => case for case in topology_study_cases(settings))
haskey(cases, case_name) || error("unknown study case: $case_name")
case = cases[case_name]

root = get(ENV, "BEAM_STUDY_OUTPUT",
           joinpath(@__DIR__, "results", "lsf_shards"))
directory = joinpath(root,
    "$(method)_$(case_name)_$(lpad(shard_index, 3, '0'))")
mkpath(directory)
record_environment(directory; settings)

points = settings["evaluation_points"]
rows = if method == :method1
    topologies = enumerate_topologies(; n=case.node_count,
        clamp_nodes=case.clamp_nodes, branch_nodes=case.branch_nodes,
        minimum_branch_degree=case.minimum_branch_degree)
    if hasproperty(case, :ignored_edges)
        topologies = filter(t -> all(!t.mask[i] for i in case.ignored_edges),
                            topologies)
    end
    selected = topologies[shard_index:shard_count:end]
    isempty(selected) && error("method1 shard contains no topologies")
    seeds = get(settings, "topology_method1_seeds", settings["seeds"])
    optimize_topologies(selected, case; seeds, points, directory)
else
    seeds = get(settings, "topology_method2_seeds", collect(1:200))
    selected = seeds[shard_index:shard_count:end]
    isempty(selected) && error("method2 shard contains no seeds")
    run_method2_initializations(case; seeds=selected,
        edge_count=length(candidate_edges(case.node_count)), directory)
end

payload = (; schema_version=1, method, case_name, shard_index, shard_count, rows)
serialize(joinpath(directory, "rows.jls"), payload)
println("Completed $(method) / $(case_name) shard $shard_index of $shard_count")
