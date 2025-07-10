module BeamStructures

    using LinearAlgebra
    using Statistics
    using DifferentialEquations
    using NonlinearSolve
    using Zygote
    using SciMLSensitivity
    using Setfield
    using Random
    using Optimisers
    using ChainRulesCore: rrule,Tangent, NoTangent, ZeroTangent, @thunk ,@non_differentiable,InplaceableThunk
    using InteractiveUtils
    import ChainRulesCore as CRC

    include("Beams.jl")
    include("Boundaries.jl")
    include("Connections.jl")
    include("Structures.jl")
    include("GroundStructures.jl")
    include("utils.jl")
    # include("Bibliography.jl")
    include("ChainRulesExt.jl")
    include("Optimizations.jl")
    
    export Connections, edge_adjacence,incidence,Adj_norm, Beam,zeros,rand

    export learningrate,changenode

    export Boundary,ExtForces, Clamp, Branch, Free, Structure,residuals!,getinitials

end