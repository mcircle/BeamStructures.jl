module BeamStructures
    # using FileIO

    # import ChainRulesCore as CRC
    # using DataFrames
    using JLD2
    using XLSX
    using LinearAlgebra
    using DifferentialEquations
    using SciMLSensitivity
    # using NLsolve
    using Setfield
    using Random
    # using ComponentArrays

    
    include("Beams.jl")
    include("Types.jl")
    include("Connections.jl")
    include("Structures.jl")
    include("Bibliography.jl")
    

    export Connections, edge_adjacence,incidence,Adj_norm, initialize, Beam

    export Boundary,ExtForces, Clamp, Branch, Free, Structure,loss!

end