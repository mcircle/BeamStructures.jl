using TOML

include("Validation.jl")
include("topology_generation.jl")
include("topology_evaluation.jl")

using .Validation
using .TopologyGeneration
using .TopologyEvaluation

settings = TOML.parsefile(joinpath(@__DIR__, "config.toml"))
adapter_path = isempty(ARGS) ? nothing : abspath(ARGS[1])
output_argument = isnothing(adapter_path) ? nothing :
                  (length(ARGS) >= 2 ? ARGS[2] : nothing)
output = !isnothing(output_argument) ? abspath(output_argument) :
         joinpath(@__DIR__, "results", "topology_" * string(time_ns()))
mkpath(output)
record_environment(output; settings)

case = nothing
if !isnothing(adapter_path)
    isfile(adapter_path) || error("adapter not found: $adapter_path")
    include(adapter_path)
    isdefined(@__MODULE__, :topology_study_case) ||
        error("adapter must define topology_study_case()")
    case = topology_study_case()
end

n = isnothing(case) ? 5 :
    (hasproperty(case, :node_count) ? case.node_count : 5)
clamps = isnothing(case) ? (1, 2, 3) :
    (hasproperty(case, :clamp_nodes) ? case.clamp_nodes : (1, 2, 3))
branches = isnothing(case) ? (4, 5) :
    (hasproperty(case, :branch_nodes) ? case.branch_nodes : (4, 5))
minimum_degree = isnothing(case) ? 2 :
    (hasproperty(case, :minimum_branch_degree) ?
     case.minimum_branch_degree : 2)

topologies = enumerate_topologies(; n, clamp_nodes=clamps,
    branch_nodes=branches, minimum_branch_degree=minimum_degree)

catalog = [(topology=t.id, elements=t.elements,
            degrees=join(t.degrees, ";")) for t in topologies]
write_rows(joinpath(output, "topology_catalog.csv"), catalog)

if isnothing(case)
    println("No adapter supplied; wrote topology catalog only.")
    println("Admissible topologies: ", length(topologies))
    println("Results: ", output)
    exit()
end

method1_seeds = get(settings, "topology_method1_seeds", settings["seeds"])
method2_seeds = get(settings, "topology_method2_seeds", collect(1:200))
points = settings["evaluation_points"]
residual_limit = get(settings, "topology_residual_limit", Inf)

runs = optimize_topologies(topologies, case; seeds=method1_seeds,
                           points, directory=output)
summary = summarize_topologies(runs; residual_limit)
write_rows(joinpath(output, "$(case.name)_topology_summary.csv"), summary)

if hasproperty(case, :method2)
    method2 = run_method2_initializations(case; seeds=method2_seeds,
        edge_count=length(candidate_edges(n)), directory=output)
    compare_method2(summary, method2; directory=output, name=case.name)
end

println("Admissible topologies: ", length(topologies))
println("Results: ", output)
