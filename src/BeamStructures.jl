module BeamStructures

using LinearAlgebra
using DifferentialEquations
using DiffEqSensitivity
using NLsolve

include("Beams.jl")
include("Connections.jl")
include("Structures.jl")

end