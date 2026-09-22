using Random

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

    positions = TopologyGeneration.random_node_positions(
        MersenneTwister(42); n=5, gridsize=100)
    @test size(positions) == (2, 5)
    @test all(x -> 0 <= x <= 100 && isinteger(x), positions)
    @test length(unique(Tuple.(eachcol(positions)))) == 5
    @test positions == TopologyGeneration.random_node_positions(
        MersenneTwister(42); n=5, gridsize=100)

    metrics = TopologyEvaluation.curve_metrics(
        [1.0 0.0 2.0; 2.0 0.0 4.0],
        [1.0 0.0 1.0; 1.0 0.0 3.0], (1.0, 1.0, 2.0))
    @test metrics.objective == 0.5
    @test length(metrics.component) == 3
    @test_throws ArgumentError TopologyEvaluation.curve_metrics(
        zeros(2, 3), zeros(2, 3), (1.0, 0.0, 1.0))

    adam = TopologyEvaluation.adam_optimize(
        x -> sum(abs2, x .- 2), [0.0, 0.0]; eta=0.1, iterations=150)
    @test adam.objective < 1e-4

    mktempdir() do directory
        selected = topologies[1:2]
        fixed_mask = copy(first(selected).mask)
        case = (
            name="fixture",
            scales=(1.0, 1.0, 1.0),
            target=points -> hcat(points, zero(points), zero(points)),
            initial=(topology, rng) -> (value=[2.0],),
            optimize=(topology, parameters) ->
                (parameters=(value=[1.0],), converged=true, residual=0.0),
            evaluate=(topology, parameters, points) ->
                hcat(parameters.value[1].*points, zero(points), zero(points)),
            method2=(rng; zero_states=false) ->
                (mask=fixed_mask, converged=true, residual=0.0))

        rows = TopologyEvaluation.optimize_topologies(selected, case;
            seeds=[1, 2], points=[-1.0, 0.0, 1.0], directory)
        @test length(rows) == 4
        @test all(row -> row.objective == 0, rows)
        @test all(row -> row.initial_objective > row.objective, rows)
        @test all(row -> row.improvement == row.initial_objective, rows)
        @test all(row -> ismissing(row.volume), rows)
        summary = TopologyEvaluation.summarize_topologies(rows;
                                                           residual_limit=1e-6)
        @test length(summary) == 2
        @test all(row -> row.success_rate == 1, summary)

        method2 = TopologyEvaluation.run_method2_initializations(case;
            seeds=1:3, edge_count=10, directory)
        @test length(method2) == 3
        @test length(unique(row.topology for row in method2)) == 1
        @test all(row -> row.initialization == "random", method2)
        zero_method2 = TopologyEvaluation.run_method2_initializations(case;
            seeds=Int[], zero_state_seeds=[7], edge_count=10, directory)
        @test length(zero_method2) == 1
        @test only(zero_method2).initialization == "zero_state"
        comparison = TopologyEvaluation.compare_method2(summary, method2;
            directory, name=case.name)
        @test all(row -> row.gap == 0, comparison)
        @test all(row -> row.frequency == 3, comparison)

        inadmissible_case = (; case..., name="inadmissible",
            admissible=mask -> false)
        rejected = TopologyEvaluation.run_method2_initializations(
            inadmissible_case; seeds=[1], edge_count=10, directory)
        @test only(rejected).status == "inadmissible"

        TopologyEvaluation.write_study_inputs(selected, (case,),
            [-1.0, 0.0, 1.0], directory)
        @test countlines(joinpath(directory, "topology_catalog.csv")) == 3
        target_path = joinpath(directory, "fixture_target.csv")
        @test countlines(target_path) == 4
        @test first(readlines(target_path)) == "point,Fx,Fy,Mz"
    end
end


@testset "Topology study adapter specification" begin
    include(joinpath(@__DIR__, "..", "validation", "topology_adapter.jl"))
    points = collect(-10.0:1.0:10.0)
    for kind in (:linear_progressive, :saddle, :valley)
        target = target_characteristic(kind, points)
        @test size(target) == (21, 3)
        @test all(iszero, target[:, 2:3])
        @test target[11, 1] == 0
    end
    @test target_characteristic(:linear_progressive, [-10.0, 10.0])[:, 1] ==
          [-10.0, 10.0]
    @test target_characteristic(:valley, [4.0, 6.0])[1, 1] < 0
    @test target_characteristic(:valley, [4.0, 6.0])[2, 1] < 0
    @test_throws ArgumentError target_characteristic(:unknown, points)

    weights = collect(0.1:0.1:1.0)
    adjacency = weighted_adjacency(weights)
    @test issymmetric(adjacency)
    @test all(iszero, diag(adjacency))
    @test adjacency[2, 1] == 0
    @test adjacency[3, 1] == weights[2]
    @test adjacency[5, 4] == weights[end]
    active_adjacency = weighted_adjacency(weights[2:end])
    @test active_adjacency == adjacency
    @test scheduled_eta(:fixed, 1, 100, 0.01) == 0.01
    @test isfinite(scheduled_eta(:cos, 1, 100, 0.01))
    @test isfinite(scheduled_eta(:inverse_sqrt, 1, 100, 0.01))
    @test scheduled_eta(:inverse_sqrt, 200, 1000, 0.005;
        parameters=200, warmups=200) ≈ 0.005

    parameters = initial_parameters(MersenneTwister(7), [-10.0f0, 0.0f0, 10.0f0])
    @test keys(parameters.nodes) == NODE_NAMES
    @test keys(parameters.beams) == BEAM_NAMES
    @test parameters.nodes.Node_1 isa BS.Clamp
    @test parameters.nodes.Node_2 isa BS.Clamp
    @test parameters.nodes.Node_3 isa BS.Branch
    @test parameters.nodes.Node_4 isa BS.Branch
    @test parameters.nodes.Node_5 isa BS.Clamp
    @test size(parameters.states) == (3, 12, 3)
    zero_parameters = initial_parameters(
        MersenneTwister(7), [-10.0f0, 0.0f0, 10.0f0]; zero_states=true)
    @test all(iszero, zero_parameters.states)

    beam = parameters.beams.Beam_1
    beam_vector = Float32[beam...]
    for fun in (BS.normfactor_m, BS.normfactor_f)
        value, pullback = CRC.rrule(fun, beam)
        @test value ≈ fun(beam)
        _, dbeam = pullback(0.7f0)
        expected = ForwardDiff.gradient(
            values -> 0.7f0 * fun(BS.Beam(values...)), beam_vector)
        @test Float32[dbeam.l, dbeam.h, dbeam.w, dbeam.E] ≈
              expected[[1, 2, 3, 5]] rtol=2f-5
    end
    norm_seed = Float32[0.2, -0.3, 0.4]
    norm_value, norm_pullback = CRC.rrule(BS.normvector, beam)
    @test norm_value ≈ BS.normvector(beam)
    _, dbeam = norm_pullback(norm_seed)
    expected = ForwardDiff.gradient(values ->
        dot(norm_seed, BS.normvector(BS.Beam(values...))), beam_vector)
    @test Float32[dbeam.l, dbeam.h, dbeam.w, dbeam.E] ≈
          expected[[1, 2, 3, 5]] rtol=2f-5

    # Exercise the same separate-argument Zygote path used by the study.
    smoke_points = Float32[0]
    smoke = initial_parameters(MersenneTwister(9), smoke_points)
    model = BS.GroundStructure()
    target = target_characteristic(:linear_progressive, smoke_points)
    value, gradients = Zygote.withgradient(
        (beams, nodes, states) -> study_loss(model, beams, nodes, states,
            ones(Float32, 10), smoke_points, target,
            (10.0f0, 10.0f0, 1000.0f0), 1.0f0),
        smoke.beams, smoke.nodes, smoke.states)
    @test isfinite(value)
    @test length(gradients) == 3
end
