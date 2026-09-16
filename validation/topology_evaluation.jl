module TopologyEvaluation

using Random, Statistics
using Optimisers, Zygote
using ..Validation: errors, write_rows
using ..TopologyGeneration: topology_id

export curve_metrics, optimize_topologies, summarize_topologies,
       compare_method2, run_method2_initializations, adam_optimize,
       write_study_inputs

"""Write the topology catalog and all target curves before optimization."""
function write_study_inputs(topologies, cases, points, directory)
    catalog = [(topology=t.id, elements=t.elements,
                degrees=join(t.degrees, ";")) for t in topologies]
    write_rows(joinpath(directory, "topology_catalog.csv"), catalog)
    for case in cases
        target = case.target(points)
        size(target) == (length(points), 3) ||
            throw(DimensionMismatch("target must have Fx, Fy, Mz columns"))
        rows = [(point=points[i], Fx=target[i, 1], Fy=target[i, 2],
                 Mz=target[i, 3]) for i in eachindex(points)]
        write_rows(joinpath(directory, "$(case.name)_target.csv"), rows)
    end
    nothing
end

"""Optimize any Optimisers-compatible parameter tree with Adam."""
function adam_optimize(loss, parameters; eta=1e-3, iterations=500,
                       callback=nothing)
    iterations >= 0 || throw(ArgumentError("iterations must be non-negative"))
    state = Optimisers.setup(Optimisers.Adam(eta), parameters)
    value = loss(parameters)
    for iteration in 1:iterations
        value, gradient = Zygote.withgradient(loss, parameters)
        isfinite(value) || error("non-finite optimization objective")
        state, parameters = Optimisers.update(state, parameters, only(gradient))
        isnothing(callback) || callback(iteration, value, parameters)
    end
    value = loss(parameters)
    (; parameters, objective=value)
end

function curve_metrics(actual, target, scales)
    size(actual) == size(target) || throw(DimensionMismatch("curve shape"))
    size(actual, 2) == 3 || throw(DimensionMismatch("expected Fx, Fy, Mz columns"))
    length(scales) == 3 && all(>(0), scales) ||
        throw(ArgumentError("three positive component scales required"))
    names = (:Fx, :Fy, :Mz)
    component = map(1:3) do j
        metric = errors(actual[:, j], target[:, j])
        (; component=names[j], metric..., normalized_mae=metric.mae/scales[j])
    end
    objective = sum(row.normalized_mae^2 for row in component)
    (; objective, component)
end

function failure_row(topology, seed, seconds, err)
    (topology=topology.id, seed, status="failed", residual=missing,
     initial_objective=missing, objective=missing, improvement=missing,
     Fx_mae=missing, Fy_mae=missing, Mz_mae=missing,
     Fx_max=missing, Fy_max=missing, Mz_max=missing, volume=missing,
     ansys_model_objective=missing, ansys_target_objective=missing,
     elements=topology.elements, seconds, message=sprint(showerror, err))
end

"""
Optimize every topology with method 1.

The case adapter supplies:
initial(topology, rng), optimize(topology, p0),
evaluate(topology, parameters, points), target(points), and scales. A volume
callback is optional. The optimizer returns
(parameters, converged, residual).
"""
function optimize_topologies(topologies, case; seeds, points, directory)
    target = case.target(points)
    size(target) == (length(points), 3) ||
        throw(DimensionMismatch("target must have Fx, Fy, Mz columns"))
    rows = NamedTuple[]
    for topology in topologies, seed in seeds
        elapsed = 0.0
        try
            initial = case.initial(topology, MersenneTwister(seed))
            before = case.evaluate(topology, initial, points)
            initial_metric = curve_metrics(before, target, case.scales)
            result = nothing
            elapsed = @elapsed result = case.optimize(topology, deepcopy(initial))
            actual = case.evaluate(topology, result.parameters, points)
            all(isfinite, actual) || error("non-finite characteristic")
            metric = curve_metrics(actual, target, case.scales)
            by_name = Dict(row.component => row for row in metric.component)
            volume = hasproperty(case, :volume) ?
                     case.volume(topology, result.parameters) : missing
            ismissing(volume) || isfinite(volume) || error("non-finite volume")
            ansys_model_objective = missing
            ansys_target_objective = missing
            if hasproperty(case, :ansys)
                ansys = case.ansys(topology, result.parameters, points)
                if !isnothing(ansys)
                    ansys_model_objective =
                        curve_metrics(ansys, actual, case.scales).objective
                    ansys_target_objective =
                        curve_metrics(ansys, target, case.scales).objective
                end
            end
            push!(rows, (topology=topology.id, seed,
                status=result.converged ? "converged" : "not_converged",
                residual=result.residual,
                initial_objective=initial_metric.objective,
                objective=metric.objective,
                improvement=initial_metric.objective-metric.objective,
                Fx_mae=by_name[:Fx].mae, Fy_mae=by_name[:Fy].mae,
                Mz_mae=by_name[:Mz].mae, Fx_max=by_name[:Fx].max_abs,
                Fy_max=by_name[:Fy].max_abs, Mz_max=by_name[:Mz].max_abs,
                volume, ansys_model_objective, ansys_target_objective,
                elements=topology.elements, seconds=elapsed, message=""))
        catch err
            push!(rows, failure_row(topology, seed, elapsed, err))
        end
    end
    write_rows(joinpath(directory, "$(case.name)_method1_runs.csv"), rows)
    rows
end

function summarize_topologies(rows; residual_limit=Inf)
    grouped = Dict{String,Vector{NamedTuple}}()
    for row in rows
        push!(get!(grouped, row.topology, NamedTuple[]), row)
    end
    summary = NamedTuple[]
    for (id, attempts) in grouped
        valid = filter(row -> row.status == "converged" &&
                              !ismissing(row.objective) &&
                              !ismissing(row.residual) &&
                              row.residual <= residual_limit, attempts)
        if isempty(valid)
            push!(summary, (topology=id, successful=0, attempts=length(attempts),
                success_rate=0.0, best_objective=missing,
                median_objective=missing, std_objective=missing,
                median_improvement=missing, best_volume=missing,
                ansys_model_objective=missing,
                ansys_target_objective=missing,
                elements=first(attempts).elements))
            continue
        end
        objectives = getproperty.(valid, :objective)
        best = valid[argmin(objectives)]
        push!(summary, (topology=id, successful=length(valid),
            attempts=length(attempts), success_rate=length(valid)/length(attempts),
            best_objective=best.objective, median_objective=median(objectives),
            std_objective=length(objectives) == 1 ? 0.0 : std(objectives),
            median_improvement=median(getproperty.(valid, :improvement)),
            best_volume=best.volume,
            ansys_model_objective=best.ansys_model_objective,
            ansys_target_objective=best.ansys_target_objective,
            elements=best.elements))
    end
    sort!(summary; by=row -> ismissing(row.best_objective) ? Inf : row.best_objective)
end

"""
Run method 2 from multiple initializations. The adapter's method2(rng) returns
(mask, converged, residual). Duplicate masks remain in the output so their
frequency and initialization sensitivity can be measured.
"""
function run_method2_initializations(case; seeds, edge_count, directory)
    rows = NamedTuple[]
    for seed in seeds
        elapsed = 0.0
        try
            result = nothing
            elapsed = @elapsed result = case.method2(MersenneTwister(seed))
            length(result.mask) == edge_count ||
                throw(DimensionMismatch("method 2 returned the wrong mask length"))
            admissible = !hasproperty(case, :admissible) ||
                         case.admissible(result.mask)
            status = !admissible ? "inadmissible" :
                     (result.converged ? "converged" : "not_converged")
            push!(rows, (seed, topology=topology_id(result.mask),
                status,
                residual=result.residual, seconds=elapsed, message=""))
        catch err
            push!(rows, (seed, topology="", status="failed", residual=missing,
                         seconds=elapsed, message=sprint(showerror, err)))
        end
    end
    write_rows(joinpath(directory, "$(case.name)_method2_runs.csv"), rows)
    rows
end

function compare_method2(summary, method2_rows; directory, name)
    lookup = Dict(row.topology => row for row in summary)
    valid_scores = [row.best_objective for row in summary
                    if !ismissing(row.best_objective)]
    best = isempty(valid_scores) ? missing : minimum(valid_scores)
    ranked = sort(unique(valid_scores))
    output = NamedTuple[]
    for row in method2_rows
        reference = get(lookup, row.topology, nothing)
        if ismissing(best) || row.status != "converged" || isnothing(reference) ||
           ismissing(reference.best_objective)
            push!(output, (seed=row.seed, topology=row.topology,
                method2_status=row.status, found_in_reference=!isnothing(reference),
                objective=missing, gap=missing, rank=missing, percentile=missing,
                frequency=count(r -> r.topology == row.topology, method2_rows)))
            continue
        end
        objective = reference.best_objective
        rank = searchsortedfirst(ranked, objective)
        percentile = length(ranked) == 1 ? 1.0 :
                     1-(rank-1)/(length(ranked)-1)
        push!(output, (seed=row.seed, topology=row.topology,
            method2_status=row.status, found_in_reference=true,
            objective, gap=objective-best, rank, percentile,
            frequency=count(r -> r.topology == row.topology, method2_rows)))
    end
    write_rows(joinpath(directory, "$(name)_method2_comparison.csv"), output)
    output
end

end
