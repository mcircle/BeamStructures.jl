using TOML

include("Validation.jl")
include("topology_generation.jl")

using .Validation: write_rows
using .TopologyGeneration: enumerate_topologies

settings = TOML.parsefile(joinpath(@__DIR__, "config.toml"))
include("topology_adapter.jl")

seed = parse(Int, get(ENV, "BEAM_DIAGNOSTIC_SEED", "11"))
iterations = parse(Int, get(ENV, "BEAM_DIAGNOSTIC_ITERATIONS", "2000"))
eta = parse(Float32, get(ENV, "BEAM_DIAGNOSTIC_ETA", "1e-2"))
log_every = parse(Int, get(ENV, "BEAM_DIAGNOSTIC_LOG_EVERY", "10"))
output = get(ENV, "BEAM_DIAGNOSTIC_OUTPUT",
             joinpath(@__DIR__, "results", "equilibrium_diagnostic.csv"))

points = Float32.(settings["evaluation_points"])
kind = Symbol(get(ENV, "BEAM_DIAGNOSTIC_CASE", "linear_progressive"))
target = Float32.(target_characteristic(kind, points;
    force_scale=settings["target_force_scale"]))
scales = Float32.((settings["target_force_scale"],
    settings["target_force_scale"],
    settings["target_force_scale"] * settings["grid_size"]))

topology = first(enumerate_topologies())
weights = Float32.(topology.mask)
initial = initial_parameters(MersenneTwister(seed), points)
model = TopologyBS.GroundStructure()
rows = NamedTuple[]

function record!(mode, iteration, state_gradient_norm, weight_gradient_norm,
                 beams, nodes, states, current_weights)
    if iteration == 0 || iteration == 1 || iteration % log_every == 0 ||
       iteration == iterations
        components = loss_components(model, beams, nodes, states,
            current_weights, points, target, scales)
        push!(rows, (; mode, iteration,
            total_loss=components.characteristic +
                Float32(settings["equilibrium_weight"]) * components.residual_mse,
            characteristic_loss=components.characteristic,
            residual_mse=components.residual_mse,
            residual_rms=components.residual_rms,
            residual_max=components.residual_max,
            state_gradient_norm,
            weight_gradient_norm))
    end
end

for mode in (:states, :beams_states, :full, :relaxed)
    parameters = deepcopy(initial)
    initial_gradient = Zygote.gradient(
        x -> equilibrium_loss(model, parameters.beams, parameters.nodes, x,
                              weights, points), parameters.states)[1]
    record!(String(mode), 0, norm(initial_gradient), NaN, parameters.beams,
            parameters.nodes, parameters.states, weights)
    callback = (iteration, _value, state_norm, weight_norm, beams, nodes,
                states, current_weights) -> record!(
                    String(mode), iteration, state_norm, weight_norm,
                    beams, nodes, states, current_weights)
    if mode === :states
        optimize_equilibrium_states(model, parameters.beams, parameters.nodes,
            parameters.states, weights, points; eta, iterations,
            callback=(iteration, value, state_norm, states) -> callback(
                iteration, value, state_norm, NaN, parameters.beams,
                parameters.nodes, states, weights))
    elseif mode === :beams_states
        optimize_equilibrium_beams_states(model, parameters.beams,
            parameters.nodes, parameters.states, weights, points;
            eta, iterations, callback)
    elseif mode === :full
        optimize_equilibrium_full(model, parameters.beams, parameters.nodes,
            parameters.states, weights, points; eta, iterations, callback)
    else
        raw_weights = zeros(Float32, length(EDGE_LIST))
        optimize_equilibrium_relaxed(model, parameters.beams, parameters.nodes,
            parameters.states, raw_weights, points; eta, iterations, callback)
    end
end

mkpath(dirname(output))
write_rows(output, rows)
println("Equilibrium diagnostic completed: ", output)
for mode in unique(row.mode for row in rows)
    selected = filter(row -> row.mode == mode, rows)
    println(mode, ": ", first(selected).residual_rms, " -> ",
            last(selected).residual_rms)
end
