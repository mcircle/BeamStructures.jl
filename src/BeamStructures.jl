module BeamStructures

    using LinearAlgebra
    using DifferentialEquations
    using DiffEqSensitivity
    using NLsolve

    include("Beams.jl")
    include("Connections.jl")
    include("Structures.jl")

    export Connections, edge_adjacence,incidence,Adj_norm

    export Boundary, Clamp, Branch, Structure,loss!

end