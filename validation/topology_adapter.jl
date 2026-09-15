import BeamStructures as TopologyBS
using LinearAlgebra, Random, Statistics
using .TopologyGeneration: candidate_edges, random_node_positions
using .TopologyEvaluation: adam_optimize

const NODE_COUNT = 5
const CLAMP_NODES = (1, 2, 5)
const BRANCH_NODES = (3, 4)
const MOVED_NODE = 5

"""
    target_characteristic(kind, displacements; force_scale=10)

Return columns `(Fx, Fy, Mz)` for one of the three study targets. The
dimensionless coordinate is `ξ = Δx/10 mm`, so `-10:10 mm` is a 20 mm span.
The third target deliberately contains a negative-stiffness interval.
"""
function target_characteristic(kind::Symbol, displacements; force_scale=10.0)
    ξ = displacements ./ 10
    fx = if kind === :linear_progressive
        force_scale .* (0.65 .* ξ .+ 0.35 .* ξ.^3)
    elseif kind === :saddle
        force_scale .* (1.5 .* ξ .- 0.5 .* ξ.^3)
    elseif kind === :valley
        force_scale .* (ξ.^3 .- 0.55 .* ξ)
    else
        throw(ArgumentError("unknown target characteristic: $kind"))
    end
    hcat(fx, zero(fx), zero(fx))
end

function weighted_adjacency(weights)
    edges = candidate_edges(NODE_COUNT)
    [i == j ? zero(eltype(weights)) :
     weights[findfirst(==((max(i, j), min(i, j))), edges)]
     for i in 1:NODE_COUNT, j in 1:NODE_COUNT]
end

function make_nodes(fixed_positions, branch_positions, displacement=0)
    positions = hcat(fixed_positions[:, 1:2], branch_positions,
                     fixed_positions[:, 5])
    positions = positions + hcat(zeros(eltype(positions), 2, 4),
                                 [displacement, zero(displacement)])
    z = zero(eltype(positions))
    TopologyBS.prepare(
        TopologyBS.Clamp(positions[1, 1], positions[2, 1], z, z, z, z),
        TopologyBS.Clamp(positions[1, 2], positions[2, 2], z, z, z, z),
        TopologyBS.Branch(positions[1, 3], positions[2, 3], z, z, z, z),
        TopologyBS.Branch(positions[1, 4], positions[2, 4], z, z, z, z),
        TopologyBS.Clamp(positions[1, 5], positions[2, 5], z, z, z, z))[2]
end

function make_beams(positions, raw_height, raw_curvature;
                    width=5.0, youngs_modulus=2.1e5)
    beams = map(enumerate(candidate_edges(NODE_COUNT))) do (index, (j, i))
        dx = positions[1, j] - positions[1, i]
        dy = positions[2, j] - positions[2, i]
        chord = hypot(dx, dy)
        chord > 0 || error("coincident nodes")
        curvature = 1.8 * tanh(raw_curvature[index]) / chord
        angle = 2 * asin(clamp(chord * curvature / 2, -0.9, 0.9))
        length = abs(curvature) < 1e-10 ? chord : angle / curvature
        start_angle = atan(dy, dx) - angle / 2
        height = 0.1 + log1p(exp(raw_height[index]))
        TopologyBS.Beam(length, height, width, curvature;
                E=youngs_modulus, θs=start_angle)
    end
    TopologyBS.prepare(beams...)[1]
end

function physical_parameters(initial, design)
    positions = hcat(initial.fixed_positions[:, 1:2],
                     design.branch_positions, initial.fixed_positions[:, 5])
    nodes = make_nodes(initial.fixed_positions, design.branch_positions)
    beams = make_beams(positions, design.raw_height, design.raw_curvature)
    (; beams, nodes)
end

function moved_nodes(nodes, displacement)
    moved = nodes.Node_5 + (; x=displacement)
    (; nodes..., Node_5=moved)
end

function moved_reaction(solutions, beams, nodes, weights)
    beam_ids = TopologyBS.getindices(NODE_COUNT)
    incident = TopologyBS.findbeamsatnode(nodes.Node_5, MOVED_NODE, beam_ids)[1]
    reaction = mapreduce(index -> weights[index] .* TopologyBS.scaleforce(
        beams[index], solutions[[1, 5, 6], 2, index]), +, incident;
        init=zeros(eltype(solutions), 3))
    # Beam state ordering is (Mz, Fx, Fy); CSV/evaluation ordering is Fx,Fy,Mz.
    reaction[[2, 3, 1]]
end

function response(initial, design, weights, displacements)
    model = TopologyBS.GroundStructure()
    parameters = physical_parameters(initial, design)
    adjacency = weighted_adjacency(weights)
    samples = map(eachindex(displacements)) do index
        nodes = moved_nodes(parameters.nodes, displacements[index])
        solutions, beams, solved_nodes = model(
            design.states[:, :, index], parameters.beams, nodes, adjacency)
        residual = zeros(eltype(solutions), size(design.states, 1),
                         size(design.states, 2))
        TopologyBS.residuals!(residual, adjacency, solutions, beams, solved_nodes)
        reaction = moved_reaction(solutions, beams, solved_nodes, weights)
        (; reaction, residual)
    end
    actual = reduce(vcat, (permutedims(sample.reaction) for sample in samples))
    residual = sqrt(mean(abs2, reduce(vcat,
        (vec(sample.residual) for sample in samples))))
    (; actual, residual)
end

function initial_design(rng, points)
    fixed_positions = random_node_positions(rng; n=NODE_COUNT, gridsize=100)
    edges = length(candidate_edges(NODE_COUNT))
    branches = copy(fixed_positions[:, collect(BRANCH_NODES)])
    design = (
        branch_positions=branches,
        raw_height=fill(log(exp(1.0) - 1), edges),
        raw_curvature=0.05 .* randn(rng, edges),
        states=0.01 .* randn(rng, 3, length(BRANCH_NODES) + edges,
                             length(points)))
    (; fixed_positions, design)
end

function study_loss(initial, design, weights, points, target, scales,
                    residual_weight, discreteness_weight=0.0)
    result = response(initial, design, weights, points)
    characteristic = mean(abs2,
        (result.actual .- target) ./ reshape(collect(scales), 1, :))
    discreteness = mean(abs2, weights .* (1 .- weights))
    characteristic + residual_weight * result.residual^2 +
        discreteness_weight * discreteness
end

function topology_study_case(kind, settings)
    points = Float64.(settings["evaluation_points"])
    target = target_characteristic(kind, points;
        force_scale=settings["target_force_scale"])
    scales = (settings["target_force_scale"],
              settings["target_force_scale"],
              settings["target_force_scale"] * settings["grid_size"])
    iterations1 = settings["adam_method1_iterations"]
    iterations2 = settings["adam_method2_iterations"]
    eta1 = settings["adam_method1_eta"]
    eta2 = settings["adam_method2_eta"]
    residual_weight = settings["equilibrium_weight"]
    threshold = settings["topology_threshold"]

    initial = (topology, rng) -> initial_design(rng, points)
    optimize = function(topology, initial_parameters)
        weights = Float64.(topology.mask)
        loss = design -> study_loss(initial_parameters, design, weights,
                                    points, target, scales, residual_weight)
        result = adam_optimize(loss, initial_parameters.design;
                               eta=eta1, iterations=iterations1)
        final = response(initial_parameters, result.parameters, weights, points)
        parameters = (; initial_parameters.fixed_positions,
                      design=result.parameters)
        (; parameters, converged=isfinite(result.objective) &&
            final.residual <= settings["topology_residual_limit"],
            residual=final.residual)
    end
    evaluate = (topology, parameters, evaluation_points) -> begin
        evaluation_points == points ||
            throw(ArgumentError("adapter evaluates the configured displacement grid"))
        response(parameters, parameters.design, Float64.(topology.mask),
                 evaluation_points).actual
    end
    method2 = function(rng)
        initial_parameters = initial_design(rng, points)
        raw_weights = randn(rng, length(candidate_edges(NODE_COUNT)))
        parameters = (; design=initial_parameters.design, raw_weights)
        loss = parameters -> begin
            weights = 1 ./ (1 .+ exp.(-parameters.raw_weights))
            study_loss(initial_parameters, parameters.design, weights, points,
                       target, scales, residual_weight,
                       settings["discreteness_weight"])
        end
        result = adam_optimize(loss, parameters; eta=eta2,
                               iterations=iterations2)
        weights = 1 ./ (1 .+ exp.(-result.parameters.raw_weights))
        final = response(initial_parameters, result.parameters.design,
                         weights, points)
        (; mask=weights .>= threshold,
           converged=isfinite(result.objective) &&
                     final.residual <= settings["topology_residual_limit"],
           residual=final.residual)
    end

    (name=String(kind), node_count=NODE_COUNT, clamp_nodes=CLAMP_NODES,
     branch_nodes=BRANCH_NODES, minimum_branch_degree=2, scales,
     target=displacements -> target_characteristic(kind, displacements;
         force_scale=settings["target_force_scale"]),
     initial, optimize, evaluate, method2)
end

topology_study_cases(settings) =
    Tuple(topology_study_case(kind, settings) for kind in
          (:linear_progressive, :saddle, :valley))
