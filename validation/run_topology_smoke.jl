using TOML

include("Validation.jl")
include("topology_generation.jl")
include("topology_evaluation.jl")

using .Validation
using .TopologyGeneration
using .TopologyEvaluation

settings = TOML.parsefile(joinpath(@__DIR__, "config.toml"))
settings["evaluation_points"] = [-10.0, 0.0, 10.0]
settings["adam_method1_iterations"] = 2
settings["adam_method2_iterations"] = 2

include("topology_adapter.jl")

output = get(ENV, "BEAM_TOPOLOGY_SMOKE_OUTPUT",
             joinpath(@__DIR__, "results", "smoke"))
mkpath(output)
record_environment(output; settings)

topologies = enumerate_topologies()[1:2]
cases = topology_study_cases(settings)
points = settings["evaluation_points"]
write_study_inputs(topologies, cases, points, output)

all_rows = NamedTuple[]
for case in cases
    method1 = optimize_topologies(topologies, case; seeds=[11], points,
                                  directory=output)
    append!(all_rows, method1)
    summary = summarize_topologies(method1; residual_limit=Inf)
    write_rows(joinpath(output, "$(case.name)_topology_summary.csv"), summary)

    method2 = run_method2_initializations(case; seeds=[1],
        edge_count=length(candidate_edges(5)), directory=output)
    append!(all_rows, method2)
    compare_method2(summary, method2; directory=output, name=case.name)
end

failed = filter(row -> row.status == "failed", all_rows)
if !isempty(failed)
    foreach(row -> println(stderr, row), failed)
    error("topology smoke study produced $(length(failed)) failed rows")
end

println("Smoke study completed: ", output)
