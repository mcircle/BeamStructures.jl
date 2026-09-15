@testset "Effective compliance and stiffness" begin
    for T in (Float32, Float64)
        tol = T === Float32 ? 2f-5 : 1e-10
        @testset "$T analytical reference values" begin
            # Fixed DOFs are deliberately non-contiguous.
            k = Matrix(Diagonal(T[9, 4, 7, 5]))
            dofs = [1, 3] => [2]
            original = copy(k)
            @test BS.effective_compliance(k, dofs) ≈ T(1 / 4)
            @test BS.effective_stiffness(k, dofs) ≈ T(4)
            @test k == original
            @test BS.effective_compliance(k, [1, 3] => [2, 4]) ≈ T(1 / 4 + 1 / 5)
            @test BS.effective_compliance(k, (1, 3) => (2,)) ≈ T(1 / 4)

            # Internal DOF relaxation gives the Schur complement 4 - 1/2.
            coupled = T[4 1; 1 2]
            @test BS.effective_stiffness(coupled, Int[] => [1]) ≈ T(3.5)
            @test BS.effective_stiffness(coupled, [2] => [1]) ≈ T(4)
            @test BS.effective_compliance(T[4 1; 2 3], Int[] => [1]) ≈ T(0.3)

            # Signed tangent stiffness is retained, rather than clamped.
            @test BS.effective_stiffness(reshape(T[-2], 1, 1), Int[] => [1]) ≈ T(-2)

            # Euler-Bernoulli cantilever tip stiffness, ordering (theta, x, y).
            # Independent analytical matrix fixture, not an assembly validation.
            L, EI, EA = T(2), T(3), T(10)
            tip = T[4*EI/L 0 -6*EI/L^2;
                    0 EA/L 0;
                    -6*EI/L^2 0 12*EI/L^3]
            @test BS.effective_compliance(tip, Int[] => [1]) ≈ L/EI rtol=tol
            @test BS.effective_compliance(tip, Int[] => [3]) ≈ L^3/(3*EI) rtol=tol
            @test BS.effective_stiffness(tip, Int[] => [2]) ≈ EA/L rtol=tol
            @test BS.effective_stiffness(tip, [1] => [3]) ≈ 12*EI/L^3 rtol=tol

            # The existing displacement-controlled force-vector API is preserved.
            @test BS.effective_movement(k, dofs) == T[0, 4, 0, 0]
        end

        @testset "$T pullbacks against ForwardDiff" begin
            symmetric = T[8 1 0 0; 1 5 0 1; 0 0 7 0; 0 1 0 3]
            nonsymmetric = T[8 1 0 0; 0 5 0 2; 0 0 7 0; 0 1 0 3]
            for k in (symmetric, nonsymmetric), loaded in ([2], [2, 4]),
                fun in (BS.effective_compliance, BS.effective_stiffness)
                dofs = [1, 3] => loaded
                expected = ForwardDiff.gradient(x -> fun(x, dofs), k)
                value, back = CRC.rrule(fun, k, dofs)
                @test value ≈ fun(k, dofs) rtol=tol
                for seed in (T(1), T(-2), T(0.25))
                    dfun, dk, ddofs = back(seed)
                    @test dfun isa CRC.NoTangent
                    @test ddofs isa CRC.NoTangent
                    @test dk ≈ seed .* expected rtol=tol atol=tol
                    @test all(iszero, dk[[1, 3], :])
                    @test all(iszero, dk[:, [1, 3]])
                end
                @test back(CRC.ZeroTangent())[2] isa CRC.AbstractZero
                @test back(CRC.@thunk(T(2)))[2] ≈ T(2) .* expected rtol=tol atol=tol
                @test Zygote.gradient(x -> fun(x, dofs), k)[1] ≈ expected rtol=tol atol=tol
            end

            # Check projection to a structured matrix tangent through its parent.
            parent = copy(nonsymmetric)
            for fun in (BS.effective_compliance, BS.effective_stiffness)
                wrapped = x -> fun(Symmetric(x, :U), [1, 3] => [2])
                @test Zygote.gradient(wrapped, parent)[1] ≈
                      ForwardDiff.gradient(wrapped, parent) rtol=tol atol=tol
            end

            # A composed optimization objective tests the outer chain rule.
            objective = k -> (BS.effective_stiffness(k, [1, 3] => [2]) - T(2))^2
            @test Zygote.gradient(objective, nonsymmetric)[1] ≈
                  ForwardDiff.gradient(objective, nonsymmetric) rtol=tol atol=tol
        end
    end

    @testset "Invalid boundary conditions and singular systems" begin
        k = Matrix{Float64}(I, 3, 3)
        invalid = (Int[] => Int[], [1] => [1], [4] => [2],
                   Int[] => [0], [1, 1] => [2], Int[] => [2, 2],
                   Int[] => [1.5])
        for fun in (BS.effective_compliance, BS.effective_stiffness)
            for dofs in invalid
                @test_throws ArgumentError fun(k, dofs)
                @test_throws ArgumentError CRC.rrule(fun, k, dofs)
            end
            @test_throws DimensionMismatch fun(ones(2, 3), Int[] => [1])
            @test_throws SingularException fun(zeros(2, 2), Int[] => [1])
        end
        # Invertible indefinite matrix with zero directional compliance.
        indefinite = [0.0 1.0; 1.0 0.0]
        @test iszero(BS.effective_compliance(indefinite, Int[] => [1]))
        @test_throws DomainError BS.effective_stiffness(indefinite, Int[] => [1])
        @test_throws DomainError CRC.rrule(BS.effective_stiffness, indefinite, Int[] => [1])
    end
end
