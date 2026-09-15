include(joinpath(@__DIR__, "..", "validation", "Validation.jl"))

@testset "Validation reports and study orchestration" begin
    @test Validation.errors([0.0], [0.0]).relative_l2 === missing
    @test Validation.errors([3.0], [1.0]).mae == 2.0
    @test_throws ArgumentError Validation.errors([NaN], [0.0])
    @test_throws DimensionMismatch Validation.errors([1.0, 2.0], [1.0])
    rows = Validation.tangent_check(q -> [2q[1]], reshape([2.0],1,1),
                                   [0.1], ones(1,1), [1e-3,1e-4])
    @test all(r -> r.max_abs < 1e-10, rows)

    mktempdir() do dir
        # Deliberately simple fixtures test bookkeeping, not mechanical validity.
        case1 = (name="fixture", initial=rng -> [2.0],
                 optimize=p -> (parameters=[1.0], converged=true, residual=0.0),
                 evaluate=(p,x) -> hcat(p[1].*x, zero(x), zero(x)),
                 target=x -> hcat(x, zero(x), zero(x)))
        rows = Validation.run_method1(case1; seeds=[1,2], points=[-1.0,0.0,1.0], directory=dir)
        @test length(rows) == 6
        @test all(r -> r.mae == 0, rows)
        curve = joinpath(dir,"fixture_seed1_curve.csv")
        comparison = Validation.compare_csv(curve, curve, joinpath(dir,"comparison.csv"))
        @test all(r -> r.max_abs == 0, comparison)

        case2 = (name="fixture", edges=2, initial=rng -> [1.0],
                 relax=p -> (parameters=p, beta=[0.2,0.8], converged=true, residual=0.0),
                 admissible=mask -> any(!iszero,mask),
                 score=(p,beta) -> (loss=sum(abs2,beta), volume=sum(beta)),
                 refine=(p,mask) -> (parameters=p, converged=true, residual=0.0))
        rows = Validation.run_method2(case2; seeds=[1], thresholds=[0.5], directory=dir)
        @test length(rows) == 7
        @test count(r -> r.stage == "enumerated", rows) == 4
        @test count(r -> r.status == "inadmissible", rows) == 1
        @test only(filter(r -> r.stage == "discrete", rows)).status == "evaluated"
        failing = merge(case1, (optimize=p -> error("solver failed"),))
        failures = Validation.run_method1(failing; seeds=[3], points=[0.0], directory=dir)
        @test only(failures).status == "failed"
    end
end
