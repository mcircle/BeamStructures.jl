include(joinpath(@__DIR__, "..", "validation", "TopologyCases.jl"))

@testset "Method 2 topology benchmarks" begin
    for T in (Float32, Float64)
        cases = TopologyCases.topology_cases(T)
        @test length(cases) == 4
        @test length(unique(case.name for case in cases)) == length(cases)
        for case in cases
            @test all(node -> node isa Union{BS.Clamp,BS.Branch}, case.nodes)
            @test issymmetric(case.adjacency)
            @test all(iszero, diag(case.adjacency))
            @test length(case.edges) == length(case.beams) == length(case.reference)
            @test 1 <= count(case.reference) < length(case.reference)
            @test TopologyCases.admissible(case, case.reference)
            @test TopologyCases.topology_quality(case, case.reference).f1 == 1
            @test TopologyCases.topology_quality(case, case.reference).hamming == 0
            @test !TopologyCases.admissible(case, falses(length(case.edges)))
            for (edge, beam) in zip(case.edges, case.beams)
                j, i = edge
                dx = case.nodes[j].x - case.nodes[i].x
                dy = case.nodes[j].y - case.nodes[i].y
                @test beam.l ≈ hypot(dx, dy)
                @test beam.θs ≈ atan(dy, dx)
            end
            for load_case in case.load_cases
                loaded = TopologyCases.apply_boundary_condition(case, load_case)
                expected = load_case.kind === :applied_load ? BS.Branch : BS.Clamp
                @test loaded[load_case.node] isa expected
                beams, nodes = TopologyCases.prepare_load_case(case, load_case)
                @test length(beams) == length(case.edges)
                @test length(nodes) == length(case.nodes)
            end
        end
    end

    force_case = only(filter(case -> case.name == "force_path",
                             TopologyCases.topology_cases(Float64)))
    loaded = TopologyCases.apply_boundary_condition(force_case,
                                                     first(force_case.load_cases))
    @test loaded[3].fx == 10
    @test loaded[3].fy == 0
    @test loaded[3].mz == 0

    quality = TopologyCases.topology_quality(force_case,
                                             trues(length(force_case.edges)))
    @test quality.recall == 1
    @test quality.precision < 1
    @test quality.hamming > 0
    @test_throws DimensionMismatch TopologyCases.topology_quality(force_case, [true])
end
