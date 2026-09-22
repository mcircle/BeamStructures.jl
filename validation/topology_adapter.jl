import BeamStructures as TopologyBS
using LinearAlgebra, Optimisers, Random, Statistics, Zygote
using .TopologyGeneration: candidate_edges, random_node_positions

const NODE_COUNT = 5
const CLAMP_NODES = (1, 2, 5)
const BRANCH_NODES = (3, 4)
const MOVED_NODE = 5
const FIXED_STIFFNESS_DOFS = [1, 2, 3, 4, 5, 6, 13, 15]
const LOADED_STIFFNESS_DOFS = [14]
const EDGE_LIST = Tuple(candidate_edges(NODE_COUNT))
const IGNORED_EDGE_IDS = (1,)
const OPTIMIZED_EDGE_IDS = Tuple(setdiff(eachindex(EDGE_LIST), IGNORED_EDGE_IDS))
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

function full_edge_weights(weights)
    if length(weights) == length(EDGE_LIST)
        allowed = map(i -> i in IGNORED_EDGE_IDS ? zero(eltype(weights)) :
                      one(eltype(weights)), eachindex(EDGE_LIST))
        return weights .* allowed
    elseif length(weights) == length(OPTIMIZED_EDGE_IDS)
        return vcat(zero(eltype(weights)), weights)
    end
    throw(DimensionMismatch("expected 9 optimized or 10 complete edge weights"))
end

function weighted_adjacency(weights)
    complete = full_edge_weights(weights)
    bounded = clamp.(complete, zero(eltype(complete)), one(eltype(complete)))
    dropdims(sum(reshape(bounded, 1, 1, :) .* EDGE_BASIS; dims=3); dims=3)
end

binary_gaussian_penalty(weights, sigma) = mean(
    TopologyBS.gaussfilter(weights, one(eltype(weights)) / 2, sigma))

function scheduled_eta(schedule::Symbol, iteration, iterations, eta;
                       phase=0.0, base=1e-5, period=iterations,
                       parameters=200, warmups=200)
    if schedule === :fixed
        eta
    elseif schedule === :inverse_sqrt
        peak_iteration = max(1, round(Int,
            warmups^(3/2) / sqrt(parameters)))
        raw = TopologyBS.learningrate(iteration, parameters, warmups)
        normalizer = TopologyBS.learningrate(
            peak_iteration, parameters, warmups)
        eta * raw / normalizer
    elseif schedule === :cos
        actual_period = period > 0 ? period : iterations
        TopologyBS.cos_learningrate(
            iteration, min(base, eta), eta, actual_period, phase)
    else
        throw(ArgumentError("unknown learning-rate schedule: $schedule"))
    end
end

function adjust_eta(state, schedule, iteration, iterations, eta, phase;
                    schedule_options...)
    value = scheduled_eta(schedule, iteration, iterations, eta;
                          phase, schedule_options...)
    Optimisers.adjust(state; eta=convert(typeof(eta), value))
end

function target_stiffness(target, points)
    n = length(points)
    n >= 2 || throw(ArgumentError("at least two points are required"))
    map(1:n) do i
        left = max(1, i - 1)
        right = min(n, i + 1)
        (target[right, 1] - target[left, 1]) /
            (points[right] - points[left])
    end
end

function effective_stiffness_curve(model, beams, nodes, states, weights, points)
    adjacency = weighted_adjacency(weights)
    map(eachindex(points)) do index
        displaced_nodes = moved_nodes(nodes, points[index])
        solutions, solved_beams, _ = model(
            states[:, :, index], beams, displaced_nodes, adjacency)
        stiffness = TopologyBS.admittance_matrix(
            solutions, adjacency, model, solved_beams)
        TopologyBS.effective_stiffness(
            stiffness, FIXED_STIFFNESS_DOFS => LOADED_STIFFNESS_DOFS)
    end
end

function topology_stiffness_fit(model, beams, nodes, states, weights, points,
                                target)
    actual = effective_stiffness_curve(
        model, beams, nodes, states, weights, points)
    desired = target_stiffness(target, points)
    scale = max(maximum(abs, desired), eps(eltype(desired)))
    mean(abs2, (actual .- desired) ./ scale)
end

function topology_stiffness_loss(model, beams, nodes, states, weights, points,
                                 target, discreteness_weight,
                                 gaussian_sigma)
    fit = topology_stiffness_fit(
        model, beams, nodes, states, weights, points, target)
    binary = binary_gaussian_penalty(
        weights, eltype(weights)(gaussian_sigma))
    fit + discreteness_weight * binary
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
    weights = full_edge_weights(weights)
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
        residual = zeros(promote_type(eltype(solutions), eltype(weights)),
                         size(states, 1), size(states, 2))
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
        residual = zeros(promote_type(eltype(states), eltype(weights)),
                         size(states, 1), size(states, 2))
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
        residual = zeros(promote_type(eltype(states), eltype(weights)),
                         size(states, 1), size(states, 2))
        residual_ = TopologyBS.residuals!(residual, adjacency, solutions,
                                          solved_beams, solved_nodes)
        mean(abs2, residual_)
    end
    mean(losses)
end

function study_loss(model, beams, nodes, states, weights, points, target,
                    scales, residual_weight, discreteness_weight=0.0,
                    gaussian_sigma=0.15)
    adjacency = weighted_adjacency(weights)
    samples = map(eachindex(points)) do index
        displaced_nodes = moved_nodes(nodes, points[index])
        solutions, solved_beams, solved_nodes = model(
            states[:, :, index], beams, displaced_nodes, adjacency)
        residual = zeros(promote_type(eltype(states), eltype(weights)),
                         size(states, 1), size(states, 2))
        residual_ = TopologyBS.residuals!(residual, adjacency, solutions,
                                          solved_beams, solved_nodes)
        reaction = moved_reaction(solutions, solved_beams, solved_nodes, weights)
        characteristic = mean(abs2,
            (reaction .- view(target, index, :)) ./ collect(scales))
        characteristic + residual_weight * mean(abs2, residual_)
    end
    discreteness = binary_gaussian_penalty(weights,
                                            eltype(weights)(gaussian_sigma))
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

function optimize_equilibrium_beams_states(model, beams, nodes, states, weights,
                                           points; eta, iterations,
                                           callback=nothing)
    beam_state = Optimisers.setup(Optimisers.Adam(eta), beams)
    value_state = Optimisers.setup(Optimisers.Adam(eta), states)
    for iteration in 1:iterations
        value, gradients = Zygote.withgradient(
            (b, x) -> equilibrium_loss(model, b, nodes, x, weights, points),
            beams, states)
        isfinite(value) || error("non-finite equilibrium objective")
        beam_state, beams = Optimisers.update(beam_state, beams, gradients[1])
        value_state, states = Optimisers.update(value_state, states, gradients[2])
        isnothing(callback) || callback(iteration, value, norm(gradients[2]),
                                        NaN, beams, nodes, states, weights)
    end
    (; beams, nodes, states, weights,
       residual_mse=equilibrium_loss(model, beams, nodes, states, weights, points))
end

function optimize_equilibrium_full(model, beams, nodes, states, weights, points;
                                   eta, iterations, callback=nothing)
    beam_state = Optimisers.setup(Optimisers.Adam(eta), beams)
    node_state = Optimisers.setup(Optimisers.Adam(eta), nodes)
    value_state = Optimisers.setup(Optimisers.Adam(eta), states)
    for iteration in 1:iterations
        value, gradients = Zygote.withgradient(
            (b, n, x) -> equilibrium_loss(model, b, n, x, weights, points),
            beams, nodes, states)
        isfinite(value) || error("non-finite equilibrium objective")
        beam_state, beams = Optimisers.update(beam_state, beams, gradients[1])
        node_state, nodes = Optimisers.update(node_state, nodes, gradients[2])
        value_state, states = Optimisers.update(value_state, states, gradients[3])
        isnothing(callback) || callback(iteration, value, norm(gradients[3]),
                                        NaN, beams, nodes, states, weights)
    end
    (; beams, nodes, states, weights,
       residual_mse=equilibrium_loss(model, beams, nodes, states, weights, points))
end

function initial_parameters(rng, points; zero_states=false)
    positions = Float32.(random_node_positions(rng; n=NODE_COUNT, gridsize=100))
    beams = make_beams(positions, rng)
    nodes = make_nodes(positions)
    shape = (3, length(BRANCH_NODES) + length(beams), length(points))
    states = zero_states ? zeros(Float32, shape) :
             0.01f0 .* randn(rng, Float32, shape)
    (; beams, nodes, states)
end

function optimize_fixed(model, parameters, weights, points, target, scales,
                        residual_weight; eta, iterations, schedule=:fixed,
                        schedule_options=NamedTuple(),
                        learning_rates=(state=eta, beam=eta, node=eta))
    beams, nodes, states = parameters.beams, parameters.nodes, parameters.states
    beam_state = Optimisers.setup(Optimisers.Adam(eta), beams)
    node_state = Optimisers.setup(Optimisers.Adam(eta), nodes)
    value_state = Optimisers.setup(Optimisers.Adam(eta), states)

    phases = (beam=-2π/3, node=-4π/3, state=0.0)
    for iteration in 1:iterations
        beam_state = adjust_eta(beam_state, schedule, iteration, iterations,
            learning_rates.beam, phases.beam; schedule_options...)
        node_state = adjust_eta(node_state, schedule, iteration, iterations,
            learning_rates.node, phases.node; schedule_options...)
        value_state = adjust_eta(value_state, schedule, iteration, iterations,
            learning_rates.state, phases.state; schedule_options...)
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
                          eta, iterations, gaussian_sigma=0.15,
                          schedule=:fixed, schedule_options=NamedTuple(),
                          learning_rates=(state=eta, beam=eta, node=eta,
                                          adjacency=eta))
    beams, nodes, states = parameters.beams, parameters.nodes, parameters.states
    beam_state = Optimisers.setup(Optimisers.Adam(eta), beams)
    node_state = Optimisers.setup(Optimisers.Adam(eta), nodes)
    value_state = Optimisers.setup(Optimisers.Adam(eta), states)
    weights = clamp.(raw_weights, zero(eltype(raw_weights)),
                     one(eltype(raw_weights)))
    weight_state = Optimisers.setup(Optimisers.Adam(eta), weights)

    phases = (beam=-π/2, node=-π, state=0.0, weight=-3π/2)
    for iteration in 1:iterations
        beam_state = adjust_eta(beam_state, schedule, iteration, iterations,
            learning_rates.beam, phases.beam; schedule_options...)
        node_state = adjust_eta(node_state, schedule, iteration, iterations,
            learning_rates.node, phases.node; schedule_options...)
        value_state = adjust_eta(value_state, schedule, iteration, iterations,
            learning_rates.state, phases.state; schedule_options...)
        weight_state = adjust_eta(weight_state, schedule, iteration, iterations,
            learning_rates.adjacency, phases.weight; schedule_options...)
        geometry_value, geometry_gradients = Zygote.withgradient(
            (w, x, y) -> study_loss(
                model, w, x, y, weights, points, target, scales,
                residual_weight),
            beams, nodes, states)
        topology_value, topology_gradient = Zygote.withgradient(
            a -> topology_stiffness_loss(
                model, beams, nodes, states, a, points, target,
                discreteness_weight, gaussian_sigma),
            weights)
        value = geometry_value + topology_value
        gradients = geometry_gradients
        weight_gradient = only(topology_gradient)
        isfinite(value) || error("non-finite optimization objective")
        beam_state, beams = Optimisers.update(beam_state, beams, gradients[1])
        node_state, nodes = Optimisers.update(node_state, nodes, gradients[2])
        value_state, states = Optimisers.update(value_state, states, gradients[3])
        weight_state, weights = Optimisers.update(
            weight_state, weights, weight_gradient)
        weights = clamp.(weights, zero(eltype(weights)), one(eltype(weights)))
    end
    objective = study_loss(model, beams, nodes, states, weights, points, target,
                           scales, residual_weight) +
                topology_stiffness_loss(
                    model, beams, nodes, states, weights, points, target,
                    discreteness_weight, gaussian_sigma)
    (; beams, nodes, states, raw_weights=weights, weights, objective)
end

function discrete_reduction_path(model, parameters, relaxed_weights,
                                 points, target, scales, residual_weight;
                                 eta, iterations, residual_limit, schedule,
                                 schedule_options, learning_rates)
    order = sortperm(relaxed_weights)
    mask = trues(length(EDGE_LIST))
    mask[collect(IGNORED_EDGE_IDS)] .= false
    current = parameters
    rows = NamedTuple[]

    function refine(mask, current, step, removed_edge, removed_weight)
        discrete_weights = Float32.(mask)
        refined = optimize_fixed(model, current, discrete_weights, points,
            target, scales, residual_weight; eta, iterations, schedule,
            schedule_options, learning_rates)
        optimized = (; beams=refined.beams, nodes=refined.nodes,
                     states=refined.states)
        actual, residual = response(model, refined.beams, refined.nodes,
                                    refined.states, discrete_weights, points)
        normalized_mae = vec(mean(abs.(actual .- target); dims=1)) ./
                         collect(scales)
        curve_objective = sum(abs2, normalized_mae)
        stiffness_error = topology_stiffness_fit(model, refined.beams,
            refined.nodes, refined.states, discrete_weights, points, target)
        (; step, removed_edge, removed_weight, mask=copy(mask),
           parameters=optimized, residual, curve_objective,
           optimization_objective=refined.objective, stiffness_error,
           converged=isfinite(refined.objective) && residual <= residual_limit)
    end

    first_step = refine(mask, current, 0, 0, missing)
    push!(rows, first_step)
    current = first_step.parameters
    step = 0
    for local_index in order
        edge = OPTIMIZED_EDGE_IDS[local_index]
        candidate = copy(mask)
        candidate[edge] = false
        TopologyGeneration.is_admissible(candidate; n=NODE_COUNT,
            clamp_nodes=CLAMP_NODES, branch_nodes=BRANCH_NODES,
            minimum_branch_degree=2) || continue
        step += 1
        result = refine(candidate, current, step, edge,
                        relaxed_weights[local_index])
        push!(rows, result)
        mask = candidate
        current = result.parameters
    end
    rank(row) = (row.converged ? 0 : 1, row.curve_objective, row.residual)
    best = rows[argmin(rank.(rows))]
    (; rows, best)
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
    schedule = Symbol(get(ENV, "BEAM_LEARNING_RATE_SCHEDULE",
                          get(settings, "learning_rate_schedule", "fixed")))
    schedule_options = (
        base=Float32(get(settings, "learning_rate_base", 1e-5)),
        period=get(settings, "learning_rate_period", 0),
        parameters=get(settings, "learning_rate_parameters", 200),
        warmups=get(settings, "learning_rate_warmups", 200))
    learning_rates = (
        state=Float32(get(settings, "learning_rate_state_peak", 5e-3)),
        beam=Float32(get(settings, "learning_rate_beam_peak", 1e-3)),
        node=Float32(get(settings, "learning_rate_node_peak", 1e-3)),
        adjacency=Float32(get(
            settings, "learning_rate_adjacency_peak", 5e-4)))

    initial = (topology, rng) -> initial_parameters(rng, points)
    optimize = function(topology, parameters)
        weights = Float32.(topology.mask)
        result = optimize_fixed(model, parameters, weights, points, target,
            scales, residual_weight; eta=Float32(settings["adam_method1_eta"]),
            iterations=settings["adam_method1_iterations"], schedule,
            schedule_options,
            learning_rates=(state=learning_rates.state,
                            beam=learning_rates.beam,
                            node=learning_rates.node))
        final = response(model, result.beams, result.nodes, result.states,
                         weights, points)
        optimized = (; beams=result.beams, nodes=result.nodes,
                     states=result.states)
        (; parameters=optimized,
           converged=isfinite(result.objective) &&
                     final[2] <= settings["topology_residual_limit"],
           residual=final[2], optimization_objective=result.objective)
    end
    evaluate = (topology, parameters, evaluation_points) -> begin
        Float32.(evaluation_points) == points ||
            throw(ArgumentError("adapter evaluates the configured displacement grid"))
        response(model, parameters.beams, parameters.nodes, parameters.states,
                 Float32.(topology.mask), points)[1]
    end
    method2 = function(rng; zero_states=false)
        parameters = initial_parameters(rng, points; zero_states)
        raw_weights = clamp.(0.5f0 .+ 0.01f0 .* randn(
            rng, Float32, length(OPTIMIZED_EDGE_IDS)), 0f0, 1f0)
        result = optimize_relaxed(model, parameters, raw_weights, points, target,
            scales, residual_weight, Float32(settings["discreteness_weight"]);
            eta=Float32(settings["adam_method2_eta"]),
            iterations=settings["adam_method2_iterations"],
            gaussian_sigma=Float32(settings["gaussian_sigma"]), schedule,
            schedule_options, learning_rates)
        relaxed_parameters = (; beams=result.beams, nodes=result.nodes,
                              states=result.states)
        relaxed_response = response(model, result.beams, result.nodes,
                                    result.states, result.weights, points)
        relaxed_stiffness_error = topology_stiffness_fit(
            model, result.beams, result.nodes, result.states, result.weights,
            points, target)
        active_bounded = clamp.(result.weights, 0f0, 1f0)
        bounded = full_edge_weights(active_bounded)
        reduction = discrete_reduction_path(model, relaxed_parameters,
            active_bounded, points, target, scales, residual_weight;
            eta=Float32(settings["adam_method1_eta"]),
            iterations=get(settings, "reduction_iterations",
                           settings["adam_method1_iterations"]),
            residual_limit=settings["topology_residual_limit"], schedule,
            schedule_options,
            learning_rates=(state=learning_rates.state,
                            beam=learning_rates.beam,
                            node=learning_rates.node))
        best = reduction.best
        discrete_weights = Float32.(best.mask)
        (; mask=best.mask, parameters=best.parameters,
           adjacency=weighted_adjacency(discrete_weights),
           continuous_adjacency=weighted_adjacency(bounded),
           weights=bounded,
           converged=best.converged,
           residual=best.residual,
           relaxed_residual=relaxed_response[2],
           optimization_objective=best.optimization_objective,
           relaxed_optimization_objective=result.objective,
           relaxed_stiffness_error,
           discrete_stiffness_error=best.stiffness_error,
           refined_stiffness_error=best.stiffness_error,
           reduction_path=reduction.rows,
           gaussian_penalty=binary_gaussian_penalty(
               active_bounded, Float32(settings["gaussian_sigma"])),
           mean_binary_distance=mean(min.(active_bounded,
                                          1f0 .- active_bounded)),
           max_binary_distance=maximum(min.(active_bounded,
                                             1f0 .- active_bounded)))
    end
    admissible = mask -> TopologyGeneration.is_admissible(mask;
        n=NODE_COUNT, clamp_nodes=CLAMP_NODES, branch_nodes=BRANCH_NODES,
        minimum_branch_degree=2)

    (name=String(kind), node_count=NODE_COUNT, clamp_nodes=CLAMP_NODES,
     branch_nodes=BRANCH_NODES, minimum_branch_degree=2, scales,
     evaluation_points=points, ignored_edges=IGNORED_EDGE_IDS,
     target=displacements -> target_characteristic(kind, displacements;
         force_scale=settings["target_force_scale"]),
     initial, optimize, evaluate, method2, admissible)
end

topology_study_cases(settings) =
    Tuple(topology_study_case(kind, settings) for kind in
          (:linear_progressive, :saddle, :valley))
