module BeamStructures
   
    # using DataFrames
    using JLD2
    using XLSX
    using LinearAlgebra
    using DifferentialEquations
    using SciMLSensitivity
    using Setfield
    using Random

    using ChainRulesCore: rrule,Tangent, NoTangent, ZeroTangent, @thunk ,@non_differentiable
    import ChainRulesCore as CRC
    # using ComponentArrays
    
    
    include("Beams.jl")
    include("Boundaries.jl")
    include("Connections.jl")
    include("Structures.jl")
    include("Bibliography.jl")
    include("utils.jl")
    include("ChainRulesExt.jl")

    export Connections, edge_adjacence,incidence,Adj_norm, initialize, Beam

    export Boundary,ExtForces, Clamp, Branch, Free, Structure,residuals!,initialize

end