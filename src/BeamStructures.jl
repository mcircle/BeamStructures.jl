module BeamStructures
    # using FileIO

    import ChainRulesCore as CRC
    using DataFrames
    using JLD2
    using XLSX
    using LinearAlgebra
    using DifferentialEquations
    using SciMLSensitivity
    using NLsolve
    using Setfield
    using Lux 
    using Random
    using NNlib
    using ComponentArrays
    using Zygote
    using Optimisers
    
    include("Beams.jl")
    include("Types.jl")
    include("Connections.jl")
    include("Structures.jl")
    include("Bibliography.jl")
    include("NeuralNets.jl")

    export Connections, edge_adjacence,incidence,Adj_norm, initialize

    export Boundary, Clamp, Branch, Structure,loss!

end