using TOML, Random, Statistics
using CSV, CairoMakie, JLD2

include("Validation.jl")
include("topology_generation.jl")
include("topology_evaluation.jl")

using .Validation
using .TopologyGeneration
using .TopologyEvaluation

settings = TOML.parsefile(joinpath(@__DIR__, "config.toml"))
include("topology_adapter.jl")

env_bool(name, default=false) = lowercase(get(ENV, name, string(default))) in
    ("1", "true", "yes", "on")

function parse_mask(value)
    text = strip(value)
    length(text) == length(EDGE_LIST) || throw(ArgumentError(
        "topology must contain $(length(EDGE_LIST)) binary entries"))
    all(character -> character in ('0', '1'), text) || throw(ArgumentError(
        "topology must be a binary string such as 0100101101"))
    BitVector(character == '1' for character in text)
end

function topology(mask)
    case_topology = only(filter(candidate -> candidate.mask == mask,
        enumerate_topologies(; n=NODE_COUNT, clamp_nodes=CLAMP_NODES,
            branch_nodes=BRANCH_NODES, minimum_branch_degree=2)))
    all(!case_topology.mask[index] for index in IGNORED_EDGE_IDS) ||
        throw(ArgumentError("ignored clamp-clamp edges must remain inactive"))
    case_topology
end

data_value(data, name, default=missing) = get(data, string(name), default)

function best_path(input, case_name, method)
    joinpath(input, "$(case_name)_$(method)_best_solution.jld2")
end

function loaded_variant(input, case_name, method)
    path = best_path(input, case_name, method)
    isfile(path) || error("best solution not found: $path")
    data = JLD2.load(path)
    mask = parse_mask(String(data_value(data, :topology)))
    raw = (beams=data_value(data, :beams), nodes=data_value(data, :nodes),
           states=data_value(data, :solution))
    parameters = extract_topology_parameters(raw, mask)
    (; method, mask, parameters,
       adjacency=weighted_adjacency(Float32.(mask)),
       continuous_adjacency=data_value(data, :continuous_adjacency,
                                        weighted_adjacency(Float32.(mask))),
       seed=Int(data_value(data, :seed, 1)),
       initialization=String(data_value(data, :initialization, "unknown")),
       stored_converged=Bool(data_value(data, :converged, false)),
       stored_residual=data_value(data, :residual),
       stored_objective=data_value(data, :objective), source=path,
       seconds=0.0)
end

function chosen_seed(input, case_name, method, fallback)
    haskey(ENV, "BEAM_SINGLE_SEED") &&
        return parse(Int, ENV["BEAM_SINGLE_SEED"])
    path = best_path(input, case_name, method)
    isfile(path) || return fallback
    Int(data_value(JLD2.load(path), :seed, fallback))
end

function chosen_topology(input, case_name)
    haskey(ENV, "BEAM_SINGLE_TOPOLOGY") &&
        return topology(parse_mask(ENV["BEAM_SINGLE_TOPOLOGY"]))
    path = best_path(input, case_name, :method1)
    isfile(path) || error(
        "set BEAM_SINGLE_TOPOLOGY or provide a method-1 best solution in $input")
    topology(parse_mask(String(data_value(JLD2.load(path), :topology))))
end

function optimize_method1(case, input, case_name)
    selected = chosen_topology(input, case_name)
    seed = chosen_seed(input, case_name, :method1, 11)
    zero_states = env_bool("BEAM_SINGLE_ZERO_STATES")
    rng = MersenneTwister(seed)
    initial = case.initial(selected, rng; zero_states)
    result = nothing
    seconds = @elapsed result = case.optimize(selected, initial)
    (; method=:method1, mask=selected.mask, parameters=result.parameters,
       adjacency=weighted_adjacency(Float32.(selected.mask)),
       continuous_adjacency=weighted_adjacency(Float32.(selected.mask)),
       seed, initialization=zero_states ? "zero_state" : "random",
       stored_converged=result.converged,
       stored_residual=result.residual,
       stored_objective=result.optimization_objective,
       source="new optimization", seconds)
end

function optimize_method2(case, input, case_name)
    seed = chosen_seed(input, case_name, :method2, 1)
    zero_states = env_bool("BEAM_SINGLE_ZERO_STATES")
    result = nothing
    seconds = @elapsed result = case.method2(
        MersenneTwister(seed); zero_states)
    (; method=:method2, mask=BitVector(result.mask),
       parameters=result.parameters, adjacency=result.adjacency,
       continuous_adjacency=result.continuous_adjacency, seed,
       initialization=zero_states ? "zero_state" : "random",
       stored_converged=result.converged,
       stored_residual=result.residual,
       stored_objective=result.optimization_objective,
       source="new optimization", seconds)
end

function evaluate_variant(case, variant, points)
    weights = Float32.(variant.mask)
    actual, residual = response(fixed_model(variant),
        variant.parameters.beams, variant.parameters.nodes,
        variant.parameters.states, weights, points)
    all(isfinite, actual) || error("non-finite characteristic")
    isfinite(residual) || error("non-finite equilibrium residual")
    target = case.target(points)
    metric = curve_metrics(actual, target, case.scales)
    (; actual, target, residual, metric,
       converged=residual <= settings["topology_residual_limit"])
end


function fixed_model(variant)
    TopologyBS.Structure(variant.adjacency)
end

function plot_variant(path, case_name, variant, evaluation, points)
    CairoMakie.activate!()
    figure = Figure(size=(1200, 900))
    geometry = Axis(figure[1:3, 1], title="$(case_name): $(variant.method)",
                    xlabel="x [mm]", ylabel="y [mm]", aspect=DataAspect(),
                    xautolimitmargin=(0.08f0, 0.08f0),
                    yautolimitmargin=(0.08f0, 0.08f0))
    nodes = variant.parameters.nodes
    for edge in findall(variant.mask)
        to_node, from_node = EDGE_LIST[edge]
        from = nodes[from_node]
        to = nodes[to_node]
        lines!(geometry, [from.x, to.x], [from.y, to.y]; linewidth=3,
               color=:steelblue)
        text!(geometry, (from.x + to.x) / 2, (from.y + to.y) / 2;
              text="B$edge", align=(:center, :bottom), fontsize=12)
    end
    node_values = collect(values(nodes))
    scatter!(geometry, getproperty.(node_values, :x),
             getproperty.(node_values, :y); color=:black, markersize=16)
    for (index, node) in enumerate(node_values)
        text!(geometry, node.x, node.y; text="N$index",
              align=(:left, :bottom), fontsize=13)
    end

    labels = (("Fx [N]", 1), ("Fy [N]", 2), ("Mz [Nm]", 3))
    for (row, (label, column)) in enumerate(labels)
        axis = Axis(figure[row, 2], xlabel=row == 3 ? "Δx [mm]" : "",
                    ylabel=label)
        lines!(axis, points, evaluation.target[:, column]; color=:black,
               linestyle=:dash, linewidth=2, label="Soll")
        lines!(axis, points, evaluation.actual[:, column]; color=:darkorange,
               linewidth=3, label="Ist")
        row == 1 && axislegend(axis; position=:lt)
    end
    Label(figure[4, 1],
        "Topologie $(topology_id(variant.mask)) | Residuum $(evaluation.residual)";
        tellwidth=false)
    save(path * ".png", figure; px_per_unit=2)
    save(path * ".pdf", figure)
end

function save_variant(output, case_name, variant, evaluation, points)
    method = String(variant.method)
    stem = joinpath(output, "$(case_name)_$(method)_single")
    plot_variant(stem, case_name, variant, evaluation, points)
    rows = [(displacement=points[index],
             target_Fx=evaluation.target[index, 1],
             actual_Fx=evaluation.actual[index, 1],
             target_Fy=evaluation.target[index, 2],
             actual_Fy=evaluation.actual[index, 2],
             target_Mz=evaluation.target[index, 3],
             actual_Mz=evaluation.actual[index, 3])
            for index in eachindex(points)]
    CSV.write(stem * "_curve.csv", rows)
    component = Dict(row.component => row for row in evaluation.metric.component)
    summary = [(case_name, method, action,
        topology=topology_id(variant.mask), elements=count(variant.mask),
        variant.seed, variant.initialization,
        residual=evaluation.residual,
        converged=evaluation.converged,
        curve_objective=evaluation.metric.objective,
        Fx_mae=component[:Fx].mae, Fy_mae=component[:Fy].mae,
        Mz_mae=component[:Mz].mae, seconds=variant.seconds,
        source=variant.source)]
    CSV.write(stem * "_summary.csv", summary)
    JLD2.jldsave(stem * "_solution.jld2";
        beams=variant.parameters.beams, nodes=variant.parameters.nodes,
        solution=variant.parameters.states, adjacency=variant.adjacency,
        continuous_adjacency=variant.continuous_adjacency,
        method, case_name, seed=variant.seed,
        initialization=variant.initialization,
        topology=topology_id(variant.mask), residual=evaluation.residual,
        objective=evaluation.metric.objective,
        converged=evaluation.converged)
    summary[1]
end

case_name = isempty(ARGS) ? "linear_progressive" : ARGS[1]
input = abspath(length(ARGS) >= 2 ? ARGS[2] :
    joinpath(@__DIR__, "results", "topology"))
output = abspath(length(ARGS) >= 3 ? ARGS[3] :
    joinpath(@__DIR__, "results", "single_$(case_name)_$(time_ns())"))
action = Symbol(lowercase(get(ENV, "BEAM_SINGLE_ACTION", "inspect")))
action in (:inspect, :optimize) ||
    error("BEAM_SINGLE_ACTION must be inspect or optimize")
method_option = lowercase(get(ENV, "BEAM_SINGLE_METHOD", "both"))
methods = method_option == "both" ? (:method1, :method2) :
          method_option == "method1" ? (:method1,) :
          method_option == "method2" ? (:method2,) :
          error("BEAM_SINGLE_METHOD must be method1, method2, or both")

cases = Dict(case.name => case for case in topology_study_cases(settings))
haskey(cases, case_name) || error(
    "unknown case $case_name; choose $(join(sort!(collect(keys(cases))), ", "))")
case = cases[case_name]
points = Float32.(settings["evaluation_points"])
mkpath(output)
record_environment(output; settings)

summaries = NamedTuple[]
for method in methods
    variant = action == :inspect ? loaded_variant(input, case_name, method) :
              method == :method1 ? optimize_method1(case, input, case_name) :
              optimize_method2(case, input, case_name)
    evaluation = evaluate_variant(case, variant, points)
    push!(summaries, save_variant(
        output, case_name, variant, evaluation, points))
end
CSV.write(joinpath(output, "single_structure_summary.csv"), summaries)

println("Single-structure $action completed: $output")
foreach(println, summaries)
