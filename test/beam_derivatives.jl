@testset "Beam RHS derivatives overwrite reusable buffers" begin
    p = BS.SciMLBase.NullParameters()
    for T in (Float32, Float64)
        state = T[0.2, 0.4, 0, 0, 0.1, -0.2, 0.3]
        reference = ForwardDiff.jacobian(x -> BS.ode(x, p, zero(T)), state)
        buffer = fill(T(NaN), 7, 7)
        BS.jac!(buffer, state, p, zero(T))
        @test buffer ≈ reference
        @test BS.jac(state, p, zero(T)) ≈ reference
        lambda = T[1, 2, 3, 4, 5, 6, 7]
        vjp = fill(T(NaN), 7)
        BS.vjp_beam!(vjp, lambda, state, zero(T))
        @test vjp ≈ -reference' * lambda
        augmented = hcat(vcat(lambda, state), vcat(2lambda, state))
        matrix_vjp = fill(T(NaN), 14, 2)
        BS.vjp!(matrix_vjp, augmented, p, zero(T))
        @test matrix_vjp[1:7,1] ≈ vjp
        @test matrix_vjp[1:7,2] ≈ 2vjp
        @test matrix_vjp[8:14,1] ≈ BS.ode(state, p, zero(T))
        @test matrix_vjp[8:14,2] ≈ matrix_vjp[8:14,1]
        curved = fill(T(NaN), 6)
        BS.vjp_curved_beam!(curved, lambda[1:6], state, zero(T))
        @test curved ≈ -reference[1:6,1:6]' * lambda[1:6]
    end
end
