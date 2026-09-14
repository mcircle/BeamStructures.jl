using Test
using LinearAlgebra
using ForwardDiff
using Zygote
import ChainRulesCore as CRC
import BeamStructures as BS

include("effective_properties.jl")

include("beam_derivatives.jl")
include("validation_environment.jl")
include("topology_cases.jl")

# Use the resolved Pkg.test environment instead of resolving/precompiling it again.
module NumericalValidationSmoke
include(joinpath(@__DIR__, "..", "validation", "run.jl"))
end
