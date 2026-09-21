using TOML

include("Validation.jl")
include("topology_generation.jl")
include("topology_evaluation.jl")

using .Validation
using .TopologyGeneration
using .TopologyEvaluation

settings = TOML.parsefile(joinpath(@__DIR__, "config.toml"))
default_adapter = joinpath(@__DIR__, "topology_adapter.jl")
mode = isempty(ARGS) ? :full :
       ARGS[1] == "--catalog-only" ? :catalog :
       ARGS[1] == "--inputs-only" ? :inputs : :full
adapter_path = mode == :catalog ? nothing :
               (mode == :full && !isempty(ARGS) ? abspath(ARGS[1]) : default_adapter)
output_argument = mode == :full ? (length(ARGS) >= 2 ? ARGS[2] : nothing) :
                  (length(ARGS) >= 2 ? ARGS[2] : nothing)
output = !isnothing(output_argument) ? abspath(output_argument) :
         joinpath(@__DIR__, "results", "topology_" * string(time_ns()))
mkpath(output)
record_environment(output; settings)

cases = ()
if !isnothing(adapter_path)
    isfile(adapter_path) || error("adapter not found: $adapter_path")
    include(adapter_path)
    isdefined(@__MODULE__, :topology_study_cases) ||
        error("adapter must define topology_study_cases(settings)")
    cases = topology_study_cases(settings)
end

n = isempty(cases) ? 5 :
    (hasproperty(first(cases), :node_count) ? first(cases).node_count : 5)
clamps = isempty(cases) ? (1, 2, 5) :
    (hasproperty(first(cases), :clamp_nodes) ? first(cases).clamp_nodes : (1, 2, 5))
branches = isempty(cases) ? (3, 4) :
    (hasproperty(first(cases), :branch_nodes) ? first(cases).branch_nodes : (3, 4))
minimum_degree = isempty(cases) ? 2 :
    (hasproperty(first(cases), :minimum_branch_degree) ?
     first(cases).minimum_branch_degree : 2)

topologies = enumerate_topologies(; n, clamp_nodes=clamps,
    branch_nodes=branches, minimum_branch_degree=minimum_degree)
if !isempty(cases) && hasproperty(first(cases), :ignored_edges)
    topologies = filter(
        t -> all(!t.mask[i] for i in first(cases).ignored_edges), topologies)
end

points = settings["evaluation_points"]
write_study_inputs(topologies, cases, points, output)

if mode != :full
    println(mode == :catalog ? "Catalog-only mode; no optimization was run." :
            "Input-only mode; no optimization was run.")
    println("Admissible topologies: ", length(topologies))
    println("Results: ", output)
    exit()
end

method1_seeds = get(settings, "topology_method1_seeds", settings["seeds"])
method2_seeds = get(settings, "topology_method2_seeds", collect(1:200))
residual_limit = get(settings, "topology_residual_limit", Inf)

for case in cases
    target = case.target(points)
    runs = optimize_topologies(topologies, case; seeds=method1_seeds,
                               points, directory=output)
    summary = summarize_topologies(runs; residual_limit)
    write_rows(joinpath(output, "$(case.name)_topology_summary.csv"), summary)

    if hasproperty(case, :method2)
        method2 = run_method2_initializations(case; seeds=method2_seeds,
            edge_count=length(candidate_edges(n)), directory=output)
        compare_method2(summary, method2; directory=output, name=case.name)
    end
end

println("Admissible topologies: ", length(topologies))
println("Results: ", output)
