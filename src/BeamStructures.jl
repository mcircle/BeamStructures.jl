module BeamStructures

    using LinearAlgebra
    using Statistics
    using DifferentialEquations
    using NonlinearSolve
    using Zygote
    using SciMLSensitivity
    using SciMLBase
    using Setfield
    using Random
    using Base
    using StaticArrays
    using Optimisers
    using ChainRulesCore: rrule,Tangent, NoTangent, ZeroTangent, @thunk ,@non_differentiable,InplaceableThunk
    using InteractiveUtils
    import ChainRulesCore as CRC

    include("BSplines.jl")
    include("Beams.jl")
    include("Boundaries.jl")
    include("Connections.jl")
    include("Structures.jl")
    include("GroundStructures.jl")
    include("utils.jl")
    # include("Bibliography.jl")
    include("Optimizations.jl")
    include("ChainRulesExt.jl")
    include("IteratorExt.jl")
    include("Exports.jl")
    
    export Connections, edge_adjacence,incidence,Adj_norm, Beam,zeros,rand,CurvedBeam,BeamElement

    export learningrate,changenode

    export Boundary,ExtForces, Clamp, Branch, Free, Structure,residuals!,getinitials,createmesh

end