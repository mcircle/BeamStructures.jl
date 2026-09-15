module Validation

using LinearAlgebra, Statistics, Random, TOML
import ForwardDiff

export errors, gradient_check, tangent_check, run_method1, run_method2,
       compare_csv, write_rows, record_environment

function errors(actual, reference)
    size(actual) == size(reference) || throw(DimensionMismatch("reference shape"))
    isempty(actual) && throw(ArgumentError("empty comparison"))
    all(isfinite, actual) && all(isfinite, reference) ||
        throw(ArgumentError("non-finite comparison"))
    difference = actual .- reference
    scale = norm(reference)
    return (mae=mean(abs, difference), max_abs=maximum(abs, difference),
            relative_l2=iszero(scale) ? missing : norm(difference)/scale)
end

function gradient_check(loss, reverse_gradient, parameters; atol, rtol)
    reference = ForwardDiff.gradient(loss, parameters)
    actual = reverse_gradient(parameters)
    e = errors(actual, reference)
    passed = maximum(abs, actual .- reference) <=
             atol + rtol * maximum(abs, reference)
    return (; e..., passed)
end

"""
Compare a tangent against central differences of independently re-equilibrated
reactions. reaction(q) must solve the SAME boundary conditions at every q.
Columns of directions carry the physical scaling of each perturbation.
"""
function tangent_check(reaction, tangent, q, directions, steps)
    rows = NamedTuple[]
    for j in axes(directions, 2), h in steps
        h > 0 || throw(ArgumentError("positive perturbation required"))
        direction = directions[:, j]
        reference = (reaction(q + h*direction) - reaction(q - h*direction))/(2*h)
        e = errors(tangent*direction, reference)
        push!(rows, (;direction=j, step=h, e...))
    end
    return rows
end

csvcell(x) = "\"" * replace(string(x), "\"" => "\"\"") * "\""
function write_rows(path, rows)
    isempty(rows) && throw(ArgumentError("no output rows"))
    mkpath(dirname(path))
    names = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(names), ","))
        for row in rows
            propertynames(row) == names || throw(ArgumentError("inconsistent columns"))
            println(io, join((csvcell(getproperty(row, n)) for n in names), ","))
        end
    end
end

function record_environment(directory; settings=Dict())
    mkpath(directory)
    root = normpath(joinpath(@__DIR__, ".."))
    commit = try
        strip(read(Cmd(["git", "-C", root, "rev-parse", "HEAD"]), String))
    catch
        "unknown"
    end
    data = Dict("julia_version"=>string(VERSION), "commit"=>commit,
                "threads"=>Threads.nthreads(), "settings"=>settings)
    open(joinpath(directory, "metadata.toml"), "w") do io
        TOML.print(io, data)
    end
    active_project = Base.active_project()
    manifest = isnothing(active_project) ? joinpath(@__DIR__, "Manifest.toml") :
               joinpath(dirname(active_project), "Manifest.toml")
    isfile(manifest) && cp(manifest, joinpath(directory, "Manifest.toml"); force=true)
end

"""
case: name, initial(rng), optimize(p0), evaluate(p, points), target(points).
optimize returns (parameters=..., converged=Bool, residual=...).
evaluate returns a matrix with columns Fx, Fy, Mz in SI units.
No solver failures are excluded from the summary.
"""
function run_method1(case; seeds, points, directory)
    rows = NamedTuple[]
    for seed in seeds
        elapsed = 0.0
        try
            initial = case.initial(MersenneTwister(seed))
            before = case.evaluate(initial, points)
            result = nothing
            elapsed = @elapsed result = case.optimize(copy(initial))
            all(isfinite, result.parameters) || error("non-finite parameters")
            after = case.evaluate(result.parameters, points)
            target = case.target(points)
            size(after) == (length(points), 3) || error("expected Fx,Fy,Mz columns")
            curve = [(point=points[i], Fx=after[i,1], Fy=after[i,2], Mz=after[i,3])
                     for i in eachindex(points)]
            write_rows(joinpath(directory, "$(case.name)_seed$(seed)_curve.csv"), curve)
            for (j, component) in enumerate(("Fx", "Fy", "Mz"))
                old = errors(before[:,j], target[:,j])
                new = errors(after[:,j], target[:,j])
                push!(rows, (seed=seed, component=component,
                    status=result.converged ? "converged" : "not_converged",
                    residual=result.residual, initial_mae=old.mae, mae=new.mae,
                    max_abs=new.max_abs, relative_l2=new.relative_l2,
                    seconds=elapsed, message=""))
            end
        catch err
            push!(rows, (seed=seed, component="all", status="failed",
                residual=missing, initial_mae=missing, mae=missing,
                max_abs=missing, relative_l2=missing, seconds=elapsed,
                message=sprint(showerror, err)))
        end
    end
    write_rows(joinpath(directory, "$(case.name)_method1.csv"), rows)
    return rows
end

"""
case: name, edges, initial(rng), relax(p0), admissible(mask),
      score(parameters, beta), refine(parameters, mask).
relax returns (parameters, beta, converged, residual).
refine returns (parameters, converged, residual).
score returns (loss, volume), using a fixed nondimensional objective.
All topology stages use identical scoring. Enumeration yields the best FOUND
reference after local parameter refinement, not a proof of global optimality.
"""
function run_method2(case; seeds, thresholds, directory, enumerate=true, max_edges=10)
    0 < case.edges <= max_edges || throw(ArgumentError("design space exceeds configured limit"))
    all(t -> 0 < t < 1, thresholds) || throw(ArgumentError("thresholds must lie in (0,1)"))
    rows = NamedTuple[]
    function record(seed, stage, threshold, mask, parameters, converged, residual, seconds)
        score = case.score(parameters, mask)
        isfinite(score.loss) && isfinite(score.volume) || error("non-finite score")
        push!(rows, (seed=seed, stage=stage, threshold=threshold,
            topology=join(mask, ";"), elements=count(>(0), mask),
            loss=score.loss, volume=score.volume, residual=residual,
            status=isnothing(converged) ? "evaluated" : (converged ? "converged" : "not_converged"), seconds=seconds, message=""))
    end
    function failure(seed, stage, threshold, mask, err; status="failed")
        push!(rows, (seed=seed, stage=stage, threshold=threshold,
            topology=join(mask, ";"), elements=count(>(0), mask),
            loss=missing, volume=missing, residual=missing, status=status,
            seconds=missing, message=string(err)))
    end
    for seed in seeds
        initial = case.initial(MersenneTwister(seed))
        try
            relaxed = nothing
            seconds = @elapsed relaxed = case.relax(copy(initial))
            length(relaxed.beta) == case.edges || error("wrong beta length")
            all(b -> isfinite(b) && 0 <= b <= 1, relaxed.beta) || error("invalid beta")
            record(seed, "relaxed", missing, relaxed.beta, relaxed.parameters,
                   relaxed.converged, relaxed.residual, seconds)
            for threshold in thresholds
                mask = Int.(relaxed.beta .>= threshold)
                if !case.admissible(mask)
                    failure(seed, "discrete", threshold, mask, "inadmissible"; status="inadmissible")
                    continue
                end
                try
                    record(seed, "discrete", threshold, mask, relaxed.parameters,
                           nothing, missing, 0.0)
                    result = nothing
                    seconds = @elapsed result = case.refine(copy(relaxed.parameters), mask)
                    record(seed, "refined", threshold, mask, result.parameters,
                           result.converged, result.residual, seconds)
                catch err
                    failure(seed, "refined", threshold, mask, sprint(showerror, err))
                end
            end
        catch err
            failure(seed, "relaxed", missing, Int[], sprint(showerror, err))
        end
        if enumerate
            for bits in 0:(2^case.edges - 1)
                mask = [Int((bits >> (i-1)) & 1) for i in 1:case.edges]
                if !case.admissible(mask)
                    failure(seed, "enumerated", missing, mask, "inadmissible"; status="inadmissible")
                    continue
                end
                try
                    result = nothing
                    seconds = @elapsed result = case.refine(copy(initial), mask)
                    record(seed, "enumerated", missing, mask, result.parameters,
                           result.converged, result.residual, seconds)
                catch err
                    failure(seed, "enumerated", missing, mask, sprint(showerror, err))
                end
            end
        end
    end
    write_rows(joinpath(directory, "$(case.name)_method2.csv"), rows)
    return rows
end

# Strict numerical CSV reader: no interpolation or silent unit conversion.
function read_curve(path)
    lines = filter(!isempty, strip.(readlines(path)))
    strip.(split(first(lines), ',')) == ["point", "Fx", "Fy", "Mz"] ||
        throw(ArgumentError("CSV header must be point,Fx,Fy,Mz"))
    data = [parse.(Float64, replace.(strip.(split(line, ',')), "\"" => ""))
            for line in lines[2:end]]
    isempty(data) && throw(ArgumentError("empty curve"))
    all(row -> length(row) == 4 && all(isfinite, row), data) || error("invalid curve")
    matrix = reduce(vcat, permutedims.(data))
    all(diff(matrix[:,1]) .> 0) || error("points must be strictly increasing")
    return matrix
end

function compare_csv(model_path, ansys_path, output; point_atol=1e-12)
    model, reference = read_curve(model_path), read_curve(ansys_path)
    size(model) == size(reference) || throw(DimensionMismatch("different sampling"))
    all(abs.(model[:,1] - reference[:,1]) .<= point_atol) ||
        error("Ansys and model must use identical sampling points")
    rows = [(component=component, errors(model[:,j+1], reference[:,j+1])...)
            for (j, component) in enumerate(("Fx", "Fy", "Mz"))]
    write_rows(output, rows)
    return rows
end

end
