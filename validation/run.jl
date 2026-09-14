using LinearAlgebra, Random, TOML
import BeamStructures as BS
import ForwardDiff, Zygote
include("Validation.jl")
using .Validation

output = get(ENV, "BEAM_VALIDATION_OUTPUT",
    length(ARGS) >= 1 ? abspath(ARGS[1]) : joinpath(@__DIR__, "results", string(time_ns())))
output = isabspath(output) ? output : joinpath(@__DIR__, "..", output)
mkpath(output)
settings = TOML.parsefile(joinpath(@__DIR__, "config.toml"))
record_environment(output; settings)
rows = NamedTuple[]
for T in (Float32, Float64)
    atol = settings[string(T)]["atol"]
    rtol = settings[string(T)]["rtol"]
    for seed in settings["seeds"]
        rng = MersenneTwister(seed)
        state = T[0.2, 0.4, 0, 0, 0.1, -0.2, 0.3] .+ T(0.05)*randn(rng, T, 7)
        p = BS.SciMLBase.NullParameters()
        reference = ForwardDiff.jacobian(x -> BS.ode(x, p, zero(T)), state)
        actual = BS.jac(state, p, zero(T))
        e = errors(actual, reference)
        push!(rows, (method="1", check="beam_rhs_jacobian", datatype=string(T), seed=seed,
                     e..., passed=isapprox(actual, reference; atol, rtol)))
        # Pre-filled storage also checks that the in-place VJP writes zero entries.
        lambda = randn(rng, T, 7)
        vjp = fill(T(NaN), 7)
        BS.vjp_beam!(vjp, lambda, state, zero(T))
        expected = -transpose(reference)*lambda
        if all(isfinite, vjp)
            e = errors(vjp, expected)
            passed = isapprox(vjp, expected; atol, rtol)
        else
            e = (mae=Inf, max_abs=Inf, relative_l2=Inf)
            passed = false
        end
        push!(rows, (method="1", check="beam_rhs_adjoint", datatype=string(T), seed=seed,
                     e..., passed=passed))

        # End-state gradients: ForwardDiff through the forward ODE versus
        # the package's backward adjoint ODE. Fixed steps avoid adaptive-path noise.
        endpoint = initial -> BS.solve(
            BS.ODEProblem(BS.ode!, initial, (zero(T), one(T)), p),
            BS.Tsit5(); adaptive=false, dt=T(0.01), save_everystep=false).u[end]
        terminal = endpoint(state)
        forward_gradient = ForwardDiff.gradient(x -> dot(lambda, endpoint(x)), state)
        backward = BS.solve(
            BS.ODEProblem(BS.vjp!, vcat(lambda,terminal), (one(T),zero(T)), p),
            BS.Tsit5(); adaptive=false, dt=-T(0.01), save_everystep=false).u[end][1:7]
        e = errors(backward, forward_gradient)
        push!(rows, (method="1", check="integrated_adjoint", datatype=string(T), seed=seed,
                     e..., passed=isapprox(backward, forward_gradient; atol, rtol)))

        A = randn(rng, T, 6, 6)
        k = transpose(A)*A + T(6)*I
        for (name, fun) in (("compliance", BS.effective_compliance),
                            ("stiffness", BS.effective_stiffness))
            result = gradient_check(x -> fun(x, [1,3] => [5]),
                x -> Zygote.gradient(y -> fun(y, [1,3] => [5]), x)[1], k; atol, rtol)
            push!(rows, (method="2", check=name*"_matrix_gradient", datatype=string(T), seed=seed,
                         result...))
        end
    end
end
write_rows(joinpath(output, "numerical_verification.csv"), rows)

# Optional project-specific adapters provide actual optimization/geometry choices.
# They run in this process and receive the common settings and output directory.
if length(ARGS) >= 2
    include(abspath(ARGS[2]))
    run_cases(settings, output)
end
println("Results: ", output)
all(row -> row.passed, rows) || error("Numerical verification failed; inspect CSV before proceeding")
