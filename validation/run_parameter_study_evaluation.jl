include("parameter_study_evaluation.jl")
using .ParameterStudyEvaluation

input = abspath(length(ARGS) >= 1 ? ARGS[1] :
    joinpath(@__DIR__, "results", "parameter_study", "aggregated"))
output = abspath(length(ARGS) >= 2 ? ARGS[2] :
    joinpath(@__DIR__, "results", "parameter_study", "evaluation"))

selected = evaluate_parameter_study(input, output)
println("Parameter-study figures and tables written to $output")
println(selected)
