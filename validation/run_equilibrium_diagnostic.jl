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
parameters = initial_parameters(MersenneTwister(seed), points)
model = TopologyBS.GroundStructure()
rows = NamedTuple[]

function record!(iteration, gradient_norm, states)
    if iteration == 0 || iteration == 1 || iteration % log_every == 0 ||
       iteration == iterations
        components = loss_components(model, parameters.beams, parameters.nodes,
            states, weights, points, target, scales)
        push!(rows, (; iteration,
            total_loss=components.characteristic +
                Float32(settings["equilibrium_weight"]) * components.residual_mse,
            characteristic_loss=components.characteristic,
            residual_mse=components.residual_mse,
            residual_rms=components.residual_rms,
            residual_max=components.residual_max,
            residual_gradient_norm=gradient_norm))
    end
end

initial_gradient = Zygote.gradient(
    x -> equilibrium_loss(model, parameters.beams, parameters.nodes, x,
                          weights, points), parameters.states)[1]
record!(0, norm(initial_gradient), parameters.states)

result = optimize_equilibrium_states(model, parameters.beams, parameters.nodes,
    parameters.states, weights, points; eta, iterations,
    callback=(iteration, _value, gradient_norm, states) ->
        record!(iteration, gradient_norm, states))

mkpath(dirname(output))
write_rows(output, rows)
println("Equilibrium diagnostic completed: ", output)
println("Initial residual RMS: ", first(rows).residual_rms)
println("Final residual RMS: ", last(rows).residual_rms)
