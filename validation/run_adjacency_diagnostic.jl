using TOML, ForwardDiff, LinearAlgebra, Random, Zygote

include("Validation.jl")
include("topology_generation.jl")

using .Validation: write_rows

settings = TOML.parsefile(joinpath(@__DIR__, "config.toml"))
include("topology_adapter.jl")

seed = parse(Int, get(ENV, "BEAM_DIAGNOSTIC_SEED", "11"))
output = get(ENV, "BEAM_ADJACENCY_DIAGNOSTIC_OUTPUT",
    joinpath(@__DIR__, "results", "adjacency_gradient_diagnostic.csv"))
points = Float32.(settings["evaluation_points"])
kind = Symbol(get(ENV, "BEAM_DIAGNOSTIC_CASE", "linear_progressive"))
target = Float32.(target_characteristic(kind, points;
    force_scale=settings["target_force_scale"]))

parameters = initial_parameters(MersenneTwister(seed), points)
model = TopologyBS.GroundStructure()
weights = fill(0.5f0, length(EDGE_LIST))
adjacency = weighted_adjacency(weights)

all(diag(adjacency) .== 0) || error("adjacency diagonal is not zero")
all(adjacency[.!Matrix{Bool}(I, NODE_COUNT, NODE_COUNT)] .== 0.5f0) ||
    error("off-diagonal adjacency entries are not 0.5")

objective(a) = topology_stiffness_loss(
    model, parameters.beams, parameters.nodes, parameters.states,
    a, points, target, 0.0f0, Float32(settings["gaussian_sigma"]))
zygote_gradient = Zygote.gradient(objective, weights)[1]
forward_gradient = ForwardDiff.gradient(objective, weights)
relative_error = norm(zygote_gradient - forward_gradient) /
                 max(norm(forward_gradient), eps(Float32))

rows = [(edge=i, node_i=EDGE_LIST[i][1], node_j=EDGE_LIST[i][2],
         zygote=zygote_gradient[i], forwarddiff=forward_gradient[i],
         absolute_error=abs(zygote_gradient[i] - forward_gradient[i]))
        for i in eachindex(weights)]
mkpath(dirname(output))
write_rows(output, rows)

println("Adjacency gradient relative error: ", relative_error)
println("Initial adjacency:\n", adjacency)
isfinite(relative_error) && relative_error <= 1.0f-4 ||
    error("adjacency stiffness gradient check failed: $relative_error")
