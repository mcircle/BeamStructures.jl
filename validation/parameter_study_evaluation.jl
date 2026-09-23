module ParameterStudyEvaluation

using CairoMakie
using CSV
using DataFrames
using LaTeXStrings
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
    save(output_path * ".pdf", figure)
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
        axislegend(axis; title="Schedule", position=:rt)
    end
    Label(figure[2, 1:3],
        "Marker: Iterationszahl (○ 500, □ 1000, △ 2000); Größe: Gesamtrechenzeit",
        tellwidth=false)
    save(output_path * ".pdf", figure)
    save(output_path * ".png", figure; px_per_unit=2)
    figure
end

function plot_pareto(pareto, output_path)
    data = finite_rows(pareto, [:elements, :objective, :residual])
    isempty(data) && return nothing
    figure = Figure(size=(760, 500))
    axis = Axis(figure[1, 1],
        xlabel=L"\mathrm{Anzahl\ der\ Balkenelemente}\;n_b",
        ylabel=L"\mathrm{Kennlinienfehler}\;J_{\mathrm{curve}}",
        title="Pareto-Lösungen der sequenziellen Topologiereduktion")
    residuals = Float64.(data.residual)
    plotted = scatter!(axis, Float64.(data.elements),
        Float64.(data.objective);
        color=residuals, colormap=:plasma, markersize=12)
    selected = data[coalesce.(data.selected, false), :]
    isempty(selected) || scatter!(axis, Float64.(selected.elements),
        Float64.(selected.objective);
        marker=:star5, markersize=22, color=:transparent,
        strokecolor=:black, strokewidth=2)
    Colorbar(figure[1, 2], plotted;
        label=L"\mathrm{Gleichgewichtsresiduum}\;r_{\mathrm{RMS}}")
    save(output_path * ".pdf", figure)
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
        score = sqrt.(
            normalized(candidates.median_objective).^2 .+
            normalized(candidates.median_residual).^2 .+
            normalized(candidates.total_seconds).^2)
        index = argmin(score)
        row = candidates[index, :]
        push!(output, (
            phase,
            config=Int(row.config),
            schedule=String(row.schedule),
            iterations=Int(row.iterations),
            learning_rate_scale=Float64(row.learning_rate_scale),
            median_objective=Float64(row.median_objective),
            median_residual=Float64(row.median_residual),
            total_seconds=Float64(row.total_seconds),
            failed=Int(row.failed),
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

    pareto_path = joinpath(input_directory, "parameter_study_pareto.csv")
    if isfile(pareto_path)
        pareto = CSV.read(pareto_path, DataFrame)
        plot_pareto(pareto, joinpath(output_directory, "reduction_pareto"))
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
