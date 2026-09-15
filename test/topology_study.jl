include(joinpath(@__DIR__, "..", "validation", "topology_generation.jl"))
include(joinpath(@__DIR__, "..", "validation", "topology_evaluation.jl"))

@testset "Topology generation and evaluation" begin
    edges = TopologyGeneration.candidate_edges(5)
    @test length(edges) == 10
    @test first(edges) == (2, 1)
    @test last(edges) == (5, 4)

    topologies = TopologyGeneration.enumerate_topologies()
    @test length(topologies) == 516
    @test all(t -> TopologyGeneration.is_admissible(t.mask), topologies)
    @test length(unique(t.id for t in topologies)) == length(topologies)
    @test all(t -> issymmetric(t.adjacency), topologies)
    @test !TopologyGeneration.is_admissible(falses(10))
    @test_throws DimensionMismatch TopologyGeneration.adjacency_matrix([true], 5)

    metrics = TopologyEvaluation.curve_metrics(
        [1.0 0.0 2.0; 2.0 0.0 4.0],
        [1.0 0.0 1.0; 1.0 0.0 3.0], (1.0, 1.0, 2.0))
    @test metrics.objective == 0.5
    @test length(metrics.component) == 3
    @test_throws ArgumentError TopologyEvaluation.curve_metrics(
        zeros(2, 3), zeros(2, 3), (1.0, 0.0, 1.0))

    mktempdir() do directory
        selected = topologies[1:2]
        fixed_mask = copy(first(selected).mask)
        case = (
            name="fixture",
            scales=(1.0, 1.0, 1.0),
            target=points -> hcat(points, zero(points), zero(points)),
            initial=(topology, rng) -> [2.0],
            optimize=(topology, parameters) ->
                (parameters=[1.0], converged=true, residual=0.0),
            evaluate=(topology, parameters, points) ->
                hcat(parameters[1].*points, zero(points), zero(points)),
            volume=(topology, parameters) -> topology.elements*parameters[1],
            method2=rng -> (mask=fixed_mask, converged=true, residual=0.0))

        rows = TopologyEvaluation.optimize_topologies(selected, case;
            seeds=[1, 2], points=[-1.0, 0.0, 1.0], directory)
        @test length(rows) == 4
        @test all(row -> row.objective == 0, rows)
        @test all(row -> row.initial_objective > row.objective, rows)
        @test all(row -> row.improvement == row.initial_objective, rows)
        summary = TopologyEvaluation.summarize_topologies(rows;
                                                           residual_limit=1e-6)
        @test length(summary) == 2
        @test all(row -> row.success_rate == 1, summary)

        method2 = TopologyEvaluation.run_method2_initializations(case;
            seeds=1:3, edge_count=10, directory)
        @test length(method2) == 3
        @test length(unique(row.topology for row in method2)) == 1
        comparison = TopologyEvaluation.compare_method2(summary, method2;
            directory, name=case.name)
        @test all(row -> row.gap == 0, comparison)
        @test all(row -> row.frequency == 3, comparison)
    end
end
