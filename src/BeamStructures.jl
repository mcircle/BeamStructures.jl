module BeamStructures
    using LinearAlgebra
    using DifferentialEquations
    using DiffEqSensitivity
    using NLsolve

    include("Beams.jl")
    include("Bibliography.jl")
    include("Connections.jl")
    include("Structures.jl")


    export Connections, edge_adjacence,incidence,Adj_norm, initialize

    export Boundary, Clamp, Branch, Structure,loss!

end