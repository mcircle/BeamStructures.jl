using Test
using LinearAlgebra
using ForwardDiff
using Zygote
import ChainRulesCore as CRC
import BeamStructures as BS

include("effective_properties.jl")

include("beam_derivatives.jl")
include("validation_environment.jl")
