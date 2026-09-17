import BeamStructures as TopologyBS
using LinearAlgebra, Optimisers, Random, Statistics, Zygote
using .TopologyGeneration: candidate_edges, random_node_positions

const NODE_COUNT = 5
const CLAMP_NODES = (1, 2, 5)
const BRANCH_NODES = (3, 4)
const MOVED_NODE = 5
const EDGE_LIST = Tuple(candidate_edges(NODE_COUNT))
const EDGE_INDEX = [i == j ? 0 :
    findfirst(==((max(i, j), min(i, j))), EDGE_LIST)
    for i in 1:NODE_COUNT, j in 1:NODE_COUNT]
const EDGE_BASIS = cat((Float32.(EDGE_INDEX .== edge)
                        for edge in eachindex(EDGE_LIST))...; dims=3)
const BEAM_NAMES = ntuple(i -> Symbol("Beam_", i), 10)
const NODE_NAMES = ntuple(i -> Symbol("Node_", i), NODE_COUNT)

"""
    target_characteristic(kind, displacements; force_scale=10)

Return columns `(Fx, Fy, Mz)` for one of the three study targets. The
dimensionless coordinate is `ξ = Δx/10 mm`, so `-10:10 mm` is a 20 mm span.
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
    dropdims(sum(reshape(weights, 1, 1, :) .* EDGE_BASIS; dims=3); dims=3)
end

"""Construct node names directly, avoiding the ambiguous internal getnames call."""
function make_nodes(positions::AbstractMatrix{T}) where {T}
    z = zero(T)
    values = (
        TopologyBS.Clamp(positions[1, 1], positions[2, 1], z, z, z, z),
        TopologyBS.Clamp(positions[1, 2], positions[2, 2], z, z, z, z),
        TopologyBS.Branch(positions[1, 3], positions[2, 3], z, z, z, z),
        TopologyBS.Branch(positions[1, 4], positions[2, 4], z, z, z, z),
        TopologyBS.Clamp(positions[1, 5], positions[2, 5], z, z, z, z))
    NamedTuple{NODE_NAMES}(values)
end

"""Construct the ten ground-structure beams without BeamStructures.prepare."""
function make_beams(positions::AbstractMatrix{T}, rng::AbstractRNG;
                    height=T(1), width=T(5), youngs_modulus=T(2.1e5)) where {T}
    values = map(EDGE_LIST) do (j, i)
        dx = positions[1, j] - positions[1, i]
        dy = positions[2, j] - positions[2, i]
        chord = hypot(dx, dy)
        chord > zero(T) || error("coincident nodes")
        curvature = T(0.05) * randn(rng, T) / chord
        angle = T(2) * asin(clamp(chord * curvature / T(2), T(-0.9), T(0.9)))
        length = abs(curvature) < sqrt(eps(T)) ? chord : angle / curvature
        start_angle = atan(dy, dx) - angle / T(2)
        TopologyBS.Beam(length, height, width, curvature;
                        E=youngs_modulus, θs=start_angle)
    end
    NamedTuple{BEAM_NAMES}(Tuple(values))
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
    # Beam state ordering is (Mz, Fx, Fy); output ordering is (Fx, Fy, Mz).
    reaction[[2, 3, 1]]
end

function response(model, beams, nodes, states, weights, displacements)
    adjacency = weighted_adjacency(weights)
    samples = map(eachindex(displacements)) do index
        displaced_nodes = moved_nodes(nodes, displacements[index])
        solutions, solved_beams, solved_nodes = model(
            states[:, :, index], beams, displaced_nodes, adjacency)
        residual = zeros(eltype(solutions), size(states, 1), size(states, 2))
        residual_ = TopologyBS.residuals!(residual, adjacency, solutions,
                                          solved_beams, solved_nodes)
        reaction = moved_reaction(solutions, solved_beams, solved_nodes, weights)
        (reaction, residual_)
    end
    actual = reduce(vcat, (permutedims(sample[1]) for sample in samples))
    residual = sqrt(mean(abs2, reduce(vcat,
        (vec(sample[2]) for sample in samples))))
    (actual, residual)
end

function loss_components(model, beams, nodes, states, weights, points, target,
                         scales)
    adjacency = weighted_adjacency(weights)
    samples = map(eachindex(points)) do index
        displaced_nodes = moved_nodes(nodes, points[index])
        solutions, solved_beams, solved_nodes = model(
            states[:, :, index], beams, displaced_nodes, adjacency)
        residual = zeros(eltype(states), size(states, 1), size(states, 2))
        residual_ = TopologyBS.residuals!(residual, adjacency, solutions,
                                          solved_beams, solved_nodes)
        reaction = moved_reaction(solutions, solved_beams, solved_nodes, weights)
        characteristic = mean(abs2,
            (reaction .- view(target, index, :)) ./ collect(scales))
        (; characteristic, residual_mse=mean(abs2, residual_),
           residual_max=maximum(abs, residual_))
    end
    (; characteristic=mean(sample.characteristic for sample in samples),
       residual_mse=mean(sample.residual_mse for sample in samples),
       residual_rms=sqrt(mean(sample.residual_mse for sample in samples)),
       residual_max=maximum(sample.residual_max for sample in samples))
end

function equilibrium_loss(model, beams, nodes, states, weights, points)
    adjacency = weighted_adjacency(weights)
    losses = map(eachindex(points)) do index
        displaced_nodes = moved_nodes(nodes, points[index])
        solutions, solved_beams, solved_nodes = model(
            states[:, :, index], beams, displaced_nodes, adjacency)
        residual = zeros(eltype(states), size(states, 1), size(states, 2))
        residual_ = TopologyBS.residuals!(residual, adjacency, solutions,
                                          solved_beams, solved_nodes)
        mean(abs2, residual_)
    end
    mean(losses)
end

function study_loss(model, beams, nodes, states, weights, points, target,
                    scales, residual_weight, discreteness_weight=0.0)
    adjacency = weighted_adjacency(weights)
    samples = map(eachindex(points)) do index
        displaced_nodes = moved_nodes(nodes, points[index])
        solutions, solved_beams, solved_nodes = model(
            states[:, :, index], beams, displaced_nodes, adjacency)
        residual = zeros(eltype(states), size(states, 1), size(states, 2))
        residual_ = TopologyBS.residuals!(residual, adjacency, solutions,
                                          solved_beams, solved_nodes)
        reaction = moved_reaction(solutions, solved_beams, solved_nodes, weights)
        characteristic = mean(abs2,
            (reaction .- view(target, index, :)) ./ collect(scales))
        characteristic + residual_weight * mean(abs2, residual_)
    end
    discreteness = mean(abs2, weights .* (one(eltype(weights)) .- weights))
    mean(samples) + discreteness_weight * discreteness
end

function optimize_equilibrium_states(model, beams, nodes, states, weights, points;
                                     eta, iterations, callback=nothing)
    optimizer_state = Optimisers.setup(Optimisers.Adam(eta), states)
    for iteration in 1:iterations
        value, gradients = Zygote.withgradient(
            x -> equilibrium_loss(model, beams, nodes, x, weights, points), states)
        isfinite(value) || error("non-finite equilibrium objective")
        gradient = only(gradients)
        gradient_norm = norm(gradient)
        optimizer_state, states = Optimisers.update(
            optimizer_state, states, gradient)
        isnothing(callback) || callback(iteration, value, gradient_norm, states)
    end
    (; states,
       residual_mse=equilibrium_loss(model, beams, nodes, states, weights, points))
end

function initial_parameters(rng, points)
    positions = Float32.(random_node_positions(rng; n=NODE_COUNT, gridsize=100))
    beams = make_beams(positions, rng)
    nodes = make_nodes(positions)
    states = 0.01f0 .* randn(rng, Float32, 3,
        length(BRANCH_NODES) + length(beams), length(points))
    (; beams, nodes, states)
end

function optimize_fixed(model, parameters, weights, points, target, scales,
                        residual_weight; eta, iterations)
    beams, nodes, states = parameters.beams, parameters.nodes, parameters.states
    beam_state = Optimisers.setup(Optimisers.Adam(eta), beams)
    node_state = Optimisers.setup(Optimisers.Adam(eta), nodes)
    value_state = Optimisers.setup(Optimisers.Adam(eta), states)

    for _ in 1:iterations
        value, gradients = Zygote.withgradient(
            (x, y, z) -> study_loss(model, x, y, z, weights, points, target,
                                     scales, residual_weight),
            beams, nodes, states)
        isfinite(value) || error("non-finite optimization objective")
        beam_state, beams = Optimisers.update(beam_state, beams, gradients[1])
        node_state, nodes = Optimisers.update(node_state, nodes, gradients[2])
        value_state, states = Optimisers.update(value_state, states, gradients[3])
    end
    objective = study_loss(model, beams, nodes, states, weights, points,
                           target, scales, residual_weight)
    (; beams, nodes, states, objective)
end

function optimize_relaxed(model, parameters, raw_weights, points, target,
                          scales, residual_weight, discreteness_weight;
                          eta, iterations)
    beams, nodes, states = parameters.beams, parameters.nodes, parameters.states
    beam_state = Optimisers.setup(Optimisers.Adam(eta), beams)
    node_state = Optimisers.setup(Optimisers.Adam(eta), nodes)
    value_state = Optimisers.setup(Optimisers.Adam(eta), states)
    weight_state = Optimisers.setup(Optimisers.Adam(eta), raw_weights)

    for _ in 1:iterations
        value, gradients = Zygote.withgradient(
            (w, x, y, z) -> begin
                weights = one(eltype(z)) ./ (one(eltype(z)) .+ exp.(-z))
                study_loss(model, w, x, y, weights, points, target, scales,
                           residual_weight, discreteness_weight)
            end,
            beams, nodes, states, raw_weights)
        isfinite(value) || error("non-finite optimization objective")
        beam_state, beams = Optimisers.update(beam_state, beams, gradients[1])
        node_state, nodes = Optimisers.update(node_state, nodes, gradients[2])
        value_state, states = Optimisers.update(value_state, states, gradients[3])
        weight_state, raw_weights = Optimisers.update(
            weight_state, raw_weights, gradients[4])
    end
    weights = one(eltype(raw_weights)) ./
              (one(eltype(raw_weights)) .+ exp.(-raw_weights))
    objective = study_loss(model, beams, nodes, states, weights, points, target,
                           scales, residual_weight, discreteness_weight)
    (; beams, nodes, states, raw_weights, weights, objective)
end

function topology_study_case(kind, settings)
    points = Float32.(settings["evaluation_points"])
    target = Float32.(target_characteristic(kind, points;
        force_scale=settings["target_force_scale"]))
    scales = Float32.((settings["target_force_scale"],
        settings["target_force_scale"],
        settings["target_force_scale"] * settings["grid_size"]))
    model = TopologyBS.GroundStructure()
    residual_weight = Float32(settings["equilibrium_weight"])

    initial = (topology, rng) -> initial_parameters(rng, points)
    optimize = function(topology, parameters)
        weights = Float32.(topology.mask)
        result = optimize_fixed(model, parameters, weights, points, target,
            scales, residual_weight; eta=Float32(settings["adam_method1_eta"]),
            iterations=settings["adam_method1_iterations"])
        final = response(model, result.beams, result.nodes, result.states,
                         weights, points)
        optimized = (; beams=result.beams, nodes=result.nodes,
                     states=result.states)
        (; parameters=optimized,
           converged=isfinite(result.objective) &&
                     final[2] <= settings["topology_residual_limit"],
           residual=final[2])
    end
    evaluate = (topology, parameters, evaluation_points) -> begin
        Float32.(evaluation_points) == points ||
            throw(ArgumentError("adapter evaluates the configured displacement grid"))
        response(model, parameters.beams, parameters.nodes, parameters.states,
                 Float32.(topology.mask), points)[1]
    end
    method2 = function(rng)
        parameters = initial_parameters(rng, points)
        raw_weights = randn(rng, Float32, length(EDGE_LIST))
        result = optimize_relaxed(model, parameters, raw_weights, points, target,
            scales, residual_weight, Float32(settings["discreteness_weight"]);
            eta=Float32(settings["adam_method2_eta"]),
            iterations=settings["adam_method2_iterations"])
        final = response(model, result.beams, result.nodes, result.states,
                         result.weights, points)
        (; mask=result.weights .>= settings["topology_threshold"],
           converged=isfinite(result.objective) &&
                     final[2] <= settings["topology_residual_limit"],
           residual=final[2])
    end
    admissible = mask -> TopologyGeneration.is_admissible(mask;
        n=NODE_COUNT, clamp_nodes=CLAMP_NODES, branch_nodes=BRANCH_NODES,
        minimum_branch_degree=2)

    (name=String(kind), node_count=NODE_COUNT, clamp_nodes=CLAMP_NODES,
     branch_nodes=BRANCH_NODES, minimum_branch_degree=2, scales,
     target=displacements -> target_characteristic(kind, displacements;
         force_scale=settings["target_force_scale"]),
     initial, optimize, evaluate, method2, admissible)
end

topology_study_cases(settings) =
    Tuple(topology_study_case(kind, settings) for kind in
          (:linear_progressive, :saddle, :valley))
