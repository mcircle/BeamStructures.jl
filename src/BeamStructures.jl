module BeamStructures
   
    # using DataFrames
    using JLD2
    using XLSX
    using LinearAlgebra
    using DifferentialEquations
    using SciMLSensitivity
    using Setfield
    using Random
    import ChainRulesCore as CRC
    # using ComponentArrays

    
    include("Beams.jl")
    include("Types.jl")
    include("Connections.jl")
    include("Structures.jl")
    include("Bibliography.jl")
    # include("ChainRulesExt.jl")

    export Connections, edge_adjacence,incidence,Adj_norm, initialize, Beam

    export Boundary,ExtForces, Clamp, Branch, Free, Structure,residuals!

end