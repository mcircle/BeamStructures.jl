include("Validation.jl")
using .Validation
length(ARGS) == 3 || error("Usage: compare_ansys.jl model.csv ansys.csv comparison.csv")
compare_csv(ARGS[1], ARGS[2], ARGS[3])
