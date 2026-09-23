module ParameterStudyEvaluation

using CairoMakie
using CSV
using DataFrames
using LaTeXStrings
using Printf
using Statistics

export evaluate_parameter_study

const PHASES = ["method1", "method2", "reduction"]
const SCHEDULES = ["fixed", "inverse_sqrt", "cos"]
const PHASE_LABELS = Dict(
    "method1" => "Methode 1",
    "method2" => "Relaxierte Methode 2",
    "reduction" => "Reduktion",
)
const SCHEDULE_COLORS = Dict(
    "fixed" => :steelblue,
    "inverse_sqrt" => :darkorange,
    "cos" => :seagreen,
)
const ITERATION_MARKERS = Dict(
    500 => :circle,
    1000 => :rect,
    2000 => :utriangle,
)

isfinitevalue(value) = !ismissing(value) && value isa Number && isfinite(value)

function require_columns(table, names)
    missing_names = setdiff(Symbol.(names), propertynames(table))
    isempty(missing_names) ||
        throw(ArgumentError("missing columns: $(join(missing_names, ", "))"))
end

function finite_rows(table, columns)
    mask = trues(nrow(table))
    for column in columns
        mask .&= map(isfinitevalue, table[!, column])
    end
    table[mask, :]
end

function log_metric(value)
    isfinitevalue(value) && value > 0 ? log10(Float64(value)) : NaN
end

function metric_range(table, metric)
    values = filter(isfinite, log_metric.(table[!, metric]))
    isempty(values) && error("no finite positive values for $metric")
    lo, hi = extrema(values)
    lo == hi ? (lo - 0.5, hi + 0.5) : (lo, hi)
end

function heatmap_values(table, phase, schedule, metric, iterations, rates)
    values = fill(NaN, length(iterations), length(rates))
    selected = table[(table.phase .== phase) .&
                     (table.schedule .== schedule), :]
    for (i, iteration) in pairs(iterations), (j, rate) in pairs(rates)
        row = selected[(selected.iterations .== iteration) .&
                       isapprox.(selected.learning_rate_scale, rate), :]
        nrow(row) == 1 || continue
        values[i, j] = log_metric(row[1, metric])
    end
    values
end

function plot_heatmaps(summary, metric, output_path; colorbar_label)
    iterations = sort(unique(Int.(summary.iterations)))
    rates = sort(unique(Float64.(summary.learning_rate_scale)))
    colorrange = metric_range(summary, metric)
    figure = Figure(size=(1180, 850))
    plotted = nothing
    for (row, phase) in pairs(PHASES), (column, schedule) in pairs(SCHEDULES)
        axis = Axis(figure[row, column],
            xlabel=L"\mathrm{Iterationszahl}\;N_{\mathrm{iter}}",
            ylabel=L"\mathrm{Lernratenfaktor}\;c_{\eta}",
            title="$(PHASE_LABELS[phase]) — $(schedule)",
            xticks=iterations, yticks=rates)
        values = heatmap_values(
            summary, phase, schedule, metric, iterations, rates)
        plotted = heatmap!(axis, iterations, rates, values;
            colormap=:viridis, colorrange)
        for i in eachindex(iterations), j in eachindex(rates)
            isfinite(values[i, j]) || continue
            text!(axis, iterations[i], rates[j];
                text=@sprintf("%.2f", values[i, j]),
                align=(:center, :center), color=:white, fontsize=11)
        end
    end
    Colorbar(figure[:, 4], plotted; label=colorbar_label)
    save(output_path * ".eps", figure)
    save(output_path * ".svg", figure)
    save(output_path * ".png", figure; px_per_unit=2)
    figure
end

function risk_heatmap_values(table, phase, schedule, metric,
                             iterations, rates, transform)
    values = fill(NaN, length(iterations), length(rates))
    selected = table[(table.phase .== phase) .&
                     (table.schedule .== schedule), :]
    for (i, iteration) in pairs(iterations), (j, rate) in pairs(rates)
        row = selected[(selected.iterations .== iteration) .&
                       isapprox.(selected.learning_rate_scale, rate), :]
        nrow(row) == 1 || continue
        value = row[1, metric]
        isfinitevalue(value) || continue
        values[i, j] = transform(Float64(value))
    end
    values
end

function plot_stiffness_risk_heatmaps(summary, metric, output_path;
                                      colorbar_label, transform=identity,
                                      annotation=value -> @sprintf("%.2f", value))
    phases = ["method2", "reduction"]
    iterations = sort(unique(Int.(summary.iterations)))
    rates = sort(unique(Float64.(summary.learning_rate_scale)))
    raw_values = [transform(Float64(value)) for value in summary[!, metric]
                  if isfinitevalue(value)]
    isempty(raw_values) && return nothing
    lo, hi = extrema(raw_values)
    colorrange = lo == hi ? (lo - 0.5, hi + 0.5) : (lo, hi)
    figure = Figure(size=(1180, 590))
    plotted = nothing
    for (row, phase) in pairs(phases), (column, schedule) in pairs(SCHEDULES)
        axis = Axis(figure[row, column],
            xlabel=L"\mathrm{Iterationszahl}\;N_{\mathrm{iter}}",
            ylabel=L"\mathrm{Lernratenfaktor}\;c_{\eta}",
            title="$(PHASE_LABELS[phase]) — $(schedule)",
            xticks=iterations, yticks=rates)
        values = risk_heatmap_values(summary, phase, schedule, metric,
                                     iterations, rates, transform)
        plotted = heatmap!(axis, iterations, rates, values;
            colormap=:magma, colorrange)
        for i in eachindex(iterations), j in eachindex(rates)
            isfinite(values[i, j]) || continue
            text!(axis, iterations[i], rates[j];
                text=annotation(values[i, j]),
                align=(:center, :center), color=:white, fontsize=11)
        end
    end
    Colorbar(figure[:, 4], plotted; label=colorbar_label)
    save(output_path * ".pdf", figure)
    save(output_path * ".svg", figure)
    save(output_path * ".eps", figure)
    save(output_path * ".png", figure; px_per_unit=2)
    figure
end

function plot_tradeoff(summary, output_path)
    data = finite_rows(summary,
        [:median_residual, :median_objective, :total_seconds])
    figure = Figure(size=(1180, 420))

    for (column, phase) in pairs(PHASES)
        axis = Axis(figure[1, column],
            xlabel=L"\mathrm{Medianes\ Residuum}\;r_{\mathrm{RMS}}",
            ylabel=L"\mathrm{Medianer\ Kennlinienfehler}\;J_{\mathrm{curve}}",
            title=PHASE_LABELS[phase], xscale=log10, yscale=log10)
        phase_data = data[data.phase .== phase, :]
        for schedule in SCHEDULES
            rows = phase_data[phase_data.schedule .== schedule, :]
            isempty(rows) && continue
            runtime = Float64.(rows.total_seconds)
            runtime_scale = maximum(runtime) == minimum(runtime) ?
                fill(14.0, length(runtime)) :
                9 .+ 13 .* (runtime .- minimum(runtime)) ./
                           (maximum(runtime) - minimum(runtime))
            scatter!(axis, Float64.(rows.median_residual),
                Float64.(rows.median_objective);
                color=SCHEDULE_COLORS[schedule],
                marker=[get(ITERATION_MARKERS, Int(value), :circle)
                        for value in rows.iterations],
                markersize=runtime_scale, label=schedule)
        end
    end
    Legend(figure[2,1],[MarkerElement(marker = :circle,color = SCHEDULE_COLORS[schedule]) for schedule in SCHEDULES],[string(s) for s in SCHEDULES], ["Schedule"],titleposition = :left,orientation = :horizontal, tellwidth=false)
    Legend(figure[2,2],[MarkerElement(marker = s,color = :transparent,strokecolor= :black,strokewidth=1) for (m,s) in ITERATION_MARKERS],[string(m) for (m,s) in ITERATION_MARKERS], ["Iterationszahlen"],orientation = :horizontal, tellwidth=false,titleposition = :left)
    Legend(figure[2,3],[MarkerElement(marker = :circle,color = :transparent,strokecolor= :black,strokewidth=1,markersize =10*mz) for mz in 1:3],["","",""], ["Gesamtrechenzeit"] ,orientation = :horizontal, tellwidth=false,titleposition = :left)
    # Label(figure[2, 3],"Größe: Gesamtrechenzeit",tellwidth=false)
    save(output_path * ".eps", figure)
    save(output_path * ".svg", figure)
    save(output_path * ".png", figure; px_per_unit=2)
    figure
end

const INITIALIZATION_LABELS = Dict(
    "random" => "Zufälliger Zustand",
    "zero_state" => "Nullzustand",
)
const INITIALIZATION_COLORS = Dict(
    "random" => :steelblue,
    "zero_state" => :darkorange,
)

function initialization_summary(runs)
    require_columns(runs, [:initialization, :objective, :residual])
    data = finite_rows(runs, [:objective, :residual])
    output = NamedTuple[]
    for group in groupby(data, :initialization)
        push!(output, (
            initialization=String(first(group.initialization)),
            runs=nrow(group),
            median_objective=median(Float64.(group.objective)),
            mean_objective=mean(Float64.(group.objective)),
            median_residual=median(Float64.(group.residual)),
            mean_residual=mean(Float64.(group.residual)),
        ))
    end
    DataFrame(output)
end

function plot_initializations(runs, output_path)
    require_columns(runs, [:initialization, :objective, :residual])
    data = finite_rows(runs, [:objective, :residual])
    isempty(data) && return nothing

    initializations = [name for name in ("random", "zero_state")
                       if name in String.(data.initialization)]
    append!(initializations,
        sort(setdiff(unique(String.(data.initialization)), initializations)))
    labels = [get(INITIALIZATION_LABELS, name, name)
              for name in initializations]

    figure = Figure(size=(900, 420))
    metrics = [
        (:objective, L"\mathrm{Kennlinienfehler}\;J_{\mathrm{curve}}"),
        (:residual, L"\mathrm{Gleichgewichtsresiduum}\;r_{\mathrm{RMS}}"),
    ]
    for (column, (metric, ylabel)) in pairs(metrics)
        axis = Axis(figure[1, column];
            xlabel="Initialisierung des Systemzustands",
            ylabel=ylabel,
            title=metric == :objective ? "Kennlinienfehler" :
                                         "Gleichgewichtsresiduum",
            yscale=log10,
            xticks=(eachindex(initializations), labels))
        for (index, initialization) in pairs(initializations)
            rows = data[String.(data.initialization) .== initialization, :]
            values = Float64.(rows[!, metric])
            offsets = 0.13 .* sin.(collect(eachindex(values)) .* 2.399963)
            color = get(INITIALIZATION_COLORS, initialization, :gray)
            scatter!(axis, index .+ offsets, values;
                color=(color, 0.28), markersize=7)
            value_median = median(values)
            lines!(axis, [index - 0.22, index + 0.22],
                [value_median, value_median];
                color=:black, linewidth=3)
            scatter!(axis, [index], [value_median];
                color=color, marker=:diamond, markersize=15,
                strokecolor=:black, strokewidth=1)
        end
    end
    Legend(figure[2, 1:2],
        [MarkerElement(marker=:circle, color=(:gray, 0.3)),
         MarkerElement(marker=:diamond, color=:gray,
                       strokecolor=:black, strokewidth=1)],
        ["einzelner Optimierungslauf", "Median"];
        orientation=:horizontal, tellwidth=false)

    save(output_path * ".pdf", figure)
    save(output_path * ".svg", figure)
    save(output_path * ".eps", figure)
    save(output_path * ".png", figure; px_per_unit=2)
    figure
end

function global_pareto_points(pareto)
    data = finite_rows(pareto, [:elements, :objective, :residual])
    isempty(data) && return data
    dominates(a, b) =
        a.elements <= b.elements &&
        a.objective <= b.objective &&
        a.residual <= b.residual &&
        (a.elements < b.elements ||
         a.objective < b.objective ||
         a.residual < b.residual)
    global_front = map(1:nrow(data)) do i
        !any(j -> j != i && dominates(data[j, :], data[i, :]),
             1:nrow(data))
    end
    data.global_pareto = global_front
    data
end

pareto_markers(data) = [
    coalesce(selected, false) ? :star5 : :circle
    for selected in data.selected
]

function plot_pareto(data, output_path)
    isempty(data) && return nothing
    figure = Figure(size=(760, 500))
    axis = Axis(figure[1, 1],
        xlabel=L"\mathrm{Anzahl\ der\ Balkenelemente}\;n_b",
        ylabel=L"\mathrm{Kennlinienfehler}\;J_{\mathrm{curve}}",
        title="Lokale und globale Pareto-Lösungen der Topologiereduktion")

    residuals = Float64.(data.residual)
    colorrange = extrema(residuals)
    colorrange = colorrange[1] == colorrange[2] ?
        (colorrange[1] - eps(), colorrange[2] + eps()) : colorrange
    local_data = data[.!data.global_pareto, :]
    global_data = data[data.global_pareto, :]

    if !isempty(local_data)
        scatter!(axis, Float64.(local_data.elements),
            Float64.(local_data.objective);
            color=Float64.(local_data.residual), colorrange,
            colormap=:plasma, marker=pareto_markers(local_data),
            markersize=[coalesce(value, false) ? 18 : 10
                        for value in local_data.selected],
            alpha=0.3)
    end
    plotted = scatter!(axis, Float64.(global_data.elements),
        Float64.(global_data.objective);
        color=Float64.(global_data.residual), colorrange,
        colormap=:plasma, marker=pareto_markers(global_data),
        markersize=[coalesce(value, false) ? 22 : 14
                    for value in global_data.selected],
        strokecolor=:black, strokewidth=1.5)

    Colorbar(figure[1, 2], plotted;
        label=L"\mathrm{Gleichgewichtsresiduum}\;r_{\mathrm{RMS}}")
    Legend(figure[2, 1:2],
        [MarkerElement(marker=:circle, color=(:gray, 0.3)),
         MarkerElement(marker=:circle, color=:gray,
                       strokecolor=:black, strokewidth=1.5),
         MarkerElement(marker=:star5, color=:gray,
                       strokecolor=:black, strokewidth=1.5)],
        ["lokal Pareto-optimal", "global Pareto-optimal",
         "ausgewählter Vorschlag"];
        orientation=:horizontal, tellwidth=false)

    save(output_path * ".pdf", figure)
    save(output_path * ".svg", figure)
    save(output_path * ".eps", figure)
    save(output_path * ".png", figure; px_per_unit=2)
    figure
end

function normalized(values)
    values = Float64.(values)
    lo, hi = extrema(values)
    hi == lo ? zeros(length(values)) : (values .- lo) ./ (hi - lo)
end

function selected_parameters(summary)
    output = NamedTuple[]
    for phase in PHASES
        candidates = finite_rows(summary[summary.phase .== phase, :],
            [:median_objective, :median_residual, :total_seconds])
        isempty(candidates) && continue
        no_failures = candidates[candidates.failed .== 0, :]
        isempty(no_failures) || (candidates = no_failures)
        score_components = [
            normalized(candidates.median_objective),
            normalized(candidates.median_residual),
            normalized(candidates.total_seconds),
        ]
        stiffness_available =
            :p90_stiffness_error in propertynames(candidates) &&
            :stiffness_error_rate_gt1 in propertynames(candidates) &&
            all(isfinitevalue, candidates.p90_stiffness_error) &&
            all(isfinitevalue, candidates.stiffness_error_rate_gt1)
        if stiffness_available
            push!(score_components,
                normalized(log1p.(Float64.(candidates.p90_stiffness_error))))
            push!(score_components,
                normalized(candidates.stiffness_error_rate_gt1))
        end
        score = sqrt.(reduce(+, component .^ 2
                             for component in score_components))
        index = argmin(score)
        row = candidates[index, :]
        push!(output, (;
            phase,
            config=Int(row.config),
            schedule=String(row.schedule),
            iterations=Int(row.iterations),
            learning_rate_scale=Float64(row.learning_rate_scale),
            median_objective=Float64(row.median_objective),
            median_residual=Float64(row.median_residual),
            total_seconds=Float64(row.total_seconds),
            failed=Int(row.failed),
            p90_stiffness_error=stiffness_available ?
                Float64(row.p90_stiffness_error) : missing,
            stiffness_error_rate_gt1=stiffness_available ?
                Float64(row.stiffness_error_rate_gt1) : missing,
            stiffness_error_rate_gt10=stiffness_available ?
                Float64(row.stiffness_error_rate_gt10) : missing,
            selection_score=score[index],
        ))
    end
    DataFrame(output)
end

function latex_escape(value)
    replace(string(value), "_" => raw"\_")
end

function write_selected_latex(path, selected)
    open(path, "w") do io
        println(io, raw"\begin{tabular}{llrrrrr}")
        println(io, raw"\toprule")
        println(io,
            raw"Phase & Schedule & $N_{\mathrm{iter}}$ & $c_\eta$ & " *
            raw"$\widetilde{J}_{\mathrm{curve}}$ & " *
            raw"$\widetilde{r}_{\mathrm{RMS}}$ & Rechenzeit [s] \\")
        println(io, raw"\midrule")
        for row in eachrow(selected)
            println(io,
                "$(latex_escape(PHASE_LABELS[row.phase])) & " *
                "$(latex_escape(row.schedule)) & $(row.iterations) & " *
                "$(round(row.learning_rate_scale; digits=2)) & " *
                "$(round(row.median_objective; sigdigits=4)) & " *
                "$(round(row.median_residual; sigdigits=4)) & " *
                "$(round(row.total_seconds; digits=1)) \\\\")
        end
        println(io, raw"\bottomrule")
        println(io, raw"\end{tabular}")
    end
end

function evaluate_parameter_study(input_directory::AbstractString,
                                  output_directory::AbstractString)
    summary_path = joinpath(input_directory, "parameter_study_summary.csv")
    isfile(summary_path) || error("missing aggregated summary: $summary_path")
    summary = CSV.read(summary_path, DataFrame)
    require_columns(summary, [
        :phase, :config, :schedule, :iterations, :learning_rate_scale,
        :median_objective, :median_residual, :total_seconds, :failed,
        :median_stiffness_error, :p90_stiffness_error,
        :maximum_stiffness_error, :stiffness_error_rate_gt1,
        :stiffness_error_rate_gt10,
    ])

    mkpath(output_directory)
    CairoMakie.activate!()
    set_theme!(theme_latexfonts())

    plot_heatmaps(summary, :median_objective,
        joinpath(output_directory, "parameter_heatmap_objective");
        colorbar_label=L"\log_{10}(\mathrm{medianer\ Kennlinienfehler})")
    plot_heatmaps(summary, :median_residual,
        joinpath(output_directory, "parameter_heatmap_residual");
        colorbar_label=L"\log_{10}(\mathrm{medianes\ Residuum})")
    plot_tradeoff(summary,
        joinpath(output_directory, "parameter_tradeoff"))
    plot_stiffness_risk_heatmaps(summary, :p90_stiffness_error,
        joinpath(output_directory, "parameter_heatmap_stiffness_p90");
        colorbar_label=L"\log_{10}(Q_{0.9}(e_k))",
        transform=log_metric)
    plot_stiffness_risk_heatmaps(summary, :stiffness_error_rate_gt1,
        joinpath(output_directory, "parameter_heatmap_stiffness_outlier_rate");
        colorbar_label=L"\mathrm{Anteil}\;e_k>1",
        annotation=value -> @sprintf("%.0f%%", 100value))

    method2_path = joinpath(input_directory, "method2_runs.csv")
    if isfile(method2_path)
        method2_runs = CSV.read(method2_path, DataFrame)
        initialization = initialization_summary(method2_runs)
        CSV.write(joinpath(output_directory, "initialization_summary.csv"),
                  initialization)
        plot_initializations(method2_runs,
            joinpath(output_directory, "initialization_comparison"))
    end

    pareto_path = joinpath(input_directory, "parameter_study_pareto.csv")
    if isfile(pareto_path)
        pareto = CSV.read(pareto_path, DataFrame)
        global_pareto = global_pareto_points(pareto)
        CSV.write(joinpath(output_directory, "global_pareto.csv"),
                  global_pareto)
        plot_pareto(global_pareto,
                    joinpath(output_directory, "reduction_pareto"))
    end

    selected = selected_parameters(summary)
    CSV.write(joinpath(output_directory, "selected_parameters.csv"), selected)
    write_selected_latex(
        joinpath(output_directory, "selected_parameters.tex"), selected)

    appendix = sort(summary,
        [:phase, :iterations, :schedule, :learning_rate_scale])
    CSV.write(joinpath(output_directory, "appendix_parameter_table.csv"),
              appendix)
    selected
end

end
