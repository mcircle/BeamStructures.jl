# Validierung von Methode 1 und 2

## Start

Mit Julia 1.12 (CI: 1.12.1, passend zum Repository-Manifest), aus dem Repository:

```sh
julia validation/setup.jl
julia --project=validation validation/run.jl
```

Die numerische Grundprüfung verwendet echte Package-Funktionen: Balken-ODE,
Jacobi-Matrix, adjungierte lokale Rückführung sowie Compliance-/Steifigkeitsgradienten.
Sie läuft für Float32/Float64 und reproduzierbare Seeds aus config.toml.
Die Matrixfälle von Methode 2 sind **numerische Referenzfälle**, keine optimierten
Balkentopologien. Ein grüner Lauf bestätigt weder die Gesamtsteifigkeitsassemblierung
noch die vollständige Adjazenz-/ODE-Ableitungskette oder die mechanische Modelltreue.

Ergebnisse stehen standardmäßig in einem neuen Unterordner von results/.
metadata.toml enthält Commit, Julia-Version, Threadzahl und Einstellungen;
das aufgelöste Manifest wird mitgesichert. Konfiguration und Einheiten vor einer
Dissertationsstudie festlegen. Toleranzen nicht allein zum Bestehen erhöhen.

## Topologiestudie auf dem 5×5-Adjazenzraum

Der eingebaute Adapter benötigt keinen Platzhalterpfad:

```sh
julia --project=validation validation/run_topology_study.jl
```

Nur den Katalog der zulässigen Topologien erzeugt:

```sh
julia --project=validation validation/run_topology_study.jl --catalog-only
```

Den Katalog und die drei Sollkennlinien ohne Optimierung als CSV erzeugen:

```sh
julia --project=validation validation/run_topology_study.jl --inputs-only validation/inputs
```

Damit entstehen `topology_catalog.csv` sowie die drei Dateien
`linear_progressive_target.csv`, `saddle_target.csv` und `valley_target.csv`.
Die erste Spalte heißt einheitlich `point`, sodass die Kennlinien direkt mit
`compare_ansys.jl` eingelesen werden können.

Ein kurzer End-to-End-Test mit einer Topologie, einem Seed, zwei
Auslenkungspunkten und je einem Adam-Schritt läuft mit:

```sh
julia --project=validation validation/run_topology_smoke.jl
```

Er verwendet dieselben Modell-, Gradienten- und CSV-Pfade wie die vollständige
Studie. Er prüft die technische Ausführbarkeit, nicht Optimierungsgüte oder
Konvergenz.

Das Gleichgewicht kann unabhängig von Geometrie und Topologie geprüft werden,
indem ausschließlich die Zustandsvariablen optimiert werden:

```sh
julia --project=validation validation/run_equilibrium_diagnostic.jl
```

Die Diagnose startet vier unabhängige Läufe von derselben Initialisierung:
nur Zustände, Zustände und Balken, Zustände/Balken/Knoten sowie zusätzlich
relaxierte Adjazenzgewichte mit dem Startwert 0,5. Die erzeugte CSV
protokolliert Gesamtverlust, Kennlinienverlust, Residuen-MSE/-RMS, maximales
absolutes Residuum sowie Zustands- und Gewichtsgradienten. Iterationszahl,
Lernrate und Ausgabe lassen sich über
`BEAM_DIAGNOSTIC_ITERATIONS`, `BEAM_DIAGNOSTIC_ETA` und
`BEAM_DIAGNOSTIC_OUTPUT` einstellen.

Der für Methode 2 relevante Adjazenzgradient wird separat geprüft:

```sh
julia --project=validation validation/run_adjacency_diagnostic.jl
```

Der Test differenziert den vorgesehenen Pfad `admittance_matrix` →
`effective_stiffness` mit Zygote und ForwardDiff. Die Sollsteifigkeit ist die
numerische Ableitung der Sollkennlinie. Die zehn unabhängigen Kanten starten
bei 0,5; daraus entsteht eine symmetrische 5×5-Adjazenzmatrix mit Nullen auf
der Diagonalen. Nach jedem Adam-Schritt werden die Kantenwerte mit `clamp` auf
`[0,1]` begrenzt. Gauß-Strafe und Binärdistanz werden als Diagnose
protokolliert, beeinflussen die Optimierung bei
`discreteness_weight = 0` jedoch nicht.

## LSF-Job-Array

Nach einmaligem Einrichten der Validierungsumgebung wird die vollständige
Studie mit folgendem Befehl eingereicht:

```sh
julia validation/setup.jl
bash validation/lsf/submit_topology_study.sh
```

Standardmäßig entstehen 48 Array-Tasks: je Sollkennlinie zwölf Shards für
Methode 1 und vier für Methode 2. Höchstens acht Tasks laufen gleichzeitig.
Die Aufteilung kann beim Einreichen angepasst werden:

```sh
METHOD1_SHARDS=16 METHOD2_SHARDS=4 MAX_CONCURRENT=8 \
  bash validation/lsf/submit_topology_study.sh
```

Jeder Task schreibt kollisionsfrei nach `validation/results/lsf_shards/`.
Nach erfolgreichem Abschluss des gesamten Arrays startet automatisch der
Merge-Job. Die finalen CSVs liegen in `validation/results/lsf_merged/`.
Für jede Sollkennlinie und Methode wird dort zusätzlich die Datei
`<kennlinie>_<methode>_best_solution.jld2` erzeugt. Sie enthält mindestens
`beams`, `nodes`, `solution` und `adjacency`; für Methode 2 wird außerdem die
kontinuierliche Adjazenzmatrix gespeichert. Nach einem vollständig erfolgreichen
Merge werden die zusammengeführten Dateien aus den eindeutig erkannten
Shard-Unterverzeichnissen standardmäßig entfernt. Mit
`BEAM_STUDY_CLEAN_SHARDS=false` bleiben die Zwischenstände erhalten. Die
LSF-Logs werden nicht automatisch gelöscht.
Schlägt ein Array-Task fehl, startet der Merge wegen der LSF-Bedingung
`done(job_id)` nicht; nach dem erneuten Ausführen fehlender Tasks kann er
manuell eingereicht werden:

```sh
bsub < validation/lsf/merge_topology_study.lsf
```

Die fünf Knoten haben fest die Rollen `Clamp, Clamp, Branch, Branch, Clamp`.
Knoten 1 und 2 sind fest; Knoten 5 wird horizontal bewegt. Für jeden Seed
werden fünf verschiedene ganzzahlige Positionen im 100×100-Raster gezogen.
Ein Seed verwendet über sämtliche festen Topologien dieselbe Startgeometrie.

Die Auslenkung läuft von -10 mm bis +10 mm. Es werden drei Sollkennlinien mit
`ξ = Δx / 10 mm` und der Kraftskala `F0` ausgewertet:

- linear-progressiv: `Fx/F0 = 0.65ξ + 0.35ξ³`
- Sattelpunkt: `Fx/F0 = 1.5ξ - 0.5ξ³`
- Tal/negative Steifigkeit: `Fx/F0 = ξ³ - 0.55ξ`

`Fy` und `Mz` sind jeweils null. Methode 1 optimiert jede graphisch zulässige
binäre Adjazenzmatrix. Methode 2 optimiert kontinuierliche Kantenwerte als
Relevanzmaß. Anschließend werden die Balken in aufsteigender Reihenfolge dieser
Werte entfernt, sofern die Struktur zulässig bleibt. Geometrie und Zustände
werden nach jedem akzeptierten Schritt erneut optimiert. Der gesamte
Reduktionspfad mit Topologie, Balkenzahl, Kennlinienfehler, Residuum und
Steifigkeitsfehler steht in `*_method2_runs.csv`. Nicht dominierte Lösungen
werden als Pareto-Lösungen markiert. Falls ein einzelner Vorschlag benötigt
wird, wird der Pareto-Punkt mit dem kleinsten normierten Abstand zum Idealpunkt
aus Kennlinienfehler, Residuum und Balkenzahl gewählt. Zusätzlich wird je
Sollkennlinie ein Lauf mit exakt null initialisierten Zuständen ausgeführt.
Beide Wege verwenden `Optimisers.Adam`; Iterationszahlen, Lernraten und
Gleichgewichtsgewicht stehen in `config.toml`. Das Volumen wird weder
berechnet noch bewertet.

Die Kante zwischen den beiden festen Einspannungen wird nicht optimiert und
bleibt in beiden Methoden null. Die adaptive Lernrate wird mit
`learning_rate_schedule` gewählt. Standardmäßig wird `inverse_sqrt` verwendet
und so normiert, dass jede Parametergruppe ihren konfigurierten Maximalwert
erreicht. Die Maximalwerte für Zustände, Balken, Knoten und Adjazenz stehen
getrennt in `config.toml`. Für einen Vergleichslauf kann die Konfiguration ohne
Dateiänderung überschrieben werden:

```sh
BEAM_LEARNING_RATE_SCHEDULE=inverse_sqrt \
BEAM_STUDY_OUTPUT=validation/results/lsf_shards_inverse \
BEAM_STUDY_MERGED_OUTPUT=validation/results/lsf_merged_inverse \
  bash validation/lsf/submit_topology_study.sh
```

Zulässige Werte sind `fixed`, `inverse_sqrt` und `cos`. Für einen belastbaren
Vergleich müssen die Varianten in getrennte Ergebnisverzeichnisse schreiben.

Im Ergebnisordner liegen `topology_catalog.csv` sowie je Kennlinie
`*_target.csv`, `*_method1_runs.csv`, `*_topology_summary.csv`,
`*_method2_runs.csv` und `*_method2_comparison.csv`. Bei einem nicht
konvergierten Referenzlauf bleiben die Vergleichswerte leer, statt den gesamten
Versuch abzubrechen.

Die Methode-2-CSV enthält zusätzlich die zehn kontinuierlichen Kantenwerte,
Gauß-Strafe, mittleren und maximalen Abstand zur Binärlösung sowie die
Steifigkeitsfehler vor Diskretisierung, direkt nach Diskretisierung und nach
der festen Nachoptimierung mit Methode 1. Nur zulässige diskrete Topologien
werden nachoptimiert und als beste Lösung berücksichtigt.

## Parameterstudie auf Batch24

Die Parameterstudie variiert getrennt für Methode 1, die relaxierte Methode 2
und die diskrete Reduktion:

- Iterationen: 500, 1000 und 2000
- Schedule: `fixed`, `inverse_sqrt` und `cos`
- Lernratenfaktor: 0,5, 1 und 2

Standardmäßig werden nur die linear-progressive Kennlinie, fünf Seeds und ein
zusätzlicher Nullzustandslauf untersucht. Die 594 logischen Teilaufgaben werden
auf einen festen Pool von 96 Array-Jobs verteilt. Jeder Worker bearbeitet
mehrere Teilaufgaben nacheinander; dadurch werden nicht hunderte einzelne Jobs
im Scheduler angelegt:

```sh
bash validation/lsf/submit_parameter_study.sh
```

Die Parallelität und der Ausgabepfad können angepasst werden:

```sh
PARAM_WORKERS=64 \
MAX_CONCURRENT=64 \
PARAM_STUDY_OUTPUT=validation/results/parameter_study_run1 \
  bash validation/lsf/submit_parameter_study.sh
```

Die Ergebnisse liegen getrennt nach Phase, Kennlinie und Parametersatz unter
`PARAM_STUDY_OUTPUT`. Die jeweils nicht untersuchten Phasen verwenden
`inverse_sqrt`, Lernratenfaktor 1 und die Basis-Iterationszahlen.

Nach Abschluss werden sämtliche Shards ohne Änderung der Rohdaten aggregiert:

```sh
bash validation/lsf/aggregate_parameter_study.sh
```

Die kompakten Dateien liegen anschließend unter
`validation/results/parameter_study/aggregated/`. Neben den drei Laufdateien
werden eine gemeinsame Parameterzusammenfassung, eine Vollständigkeitsprüfung,
die Pareto-Lösungen der Reduktionsphase und jeweils die beste JLD2-Lösung pro
Phase erzeugt. Bei fehlenden Shards endet das Skript mit Exit-Code 2 und listet
sie in `parameter_study_completeness.csv`.

Die aggregierten Daten werden mit CairoMakie und dem LaTeX-Schriftthema
ausgewertet:

```sh
julia --project=validation validation/run_parameter_study_evaluation.jl
```

Unter `validation/results/parameter_study/evaluation/` entstehen:

- Heatmaps für Kennlinienfehler und Gleichgewichtsresiduum als PDF und PNG,
- das Trade-off-Diagramm aus Kennlinienfehler und Residuum,
- das Pareto-Diagramm der Topologiereduktion,
- `selected_parameters.csv` und `selected_parameters.tex` für den Haupttext,
- `appendix_parameter_table.csv` mit allen Parameterkonfigurationen.

Alle Diagramme besitzen beschriftete Achsen und verwenden logarithmische
Darstellungen nur dort, wo dies explizit in Achse oder Farbskala angegeben ist.
Der automatisch ausgewählte Parametersatz minimiert den gleich gewichteten,
normierten Abstand aus medianem Kennlinienfehler, medianem Residuum und
Gesamtrechenzeit. Diese Auswahl ist als Vorschlag zu prüfen, nicht als
physikalisch zwingende Gewichtung.

## Eigene Optimierungsfälle

```sh
julia --project=validation validation/run.jl validation/results/study path/to/cases.jl
```

cases.jl definiert `run_cases(settings, output)` und ruft die Funktionen
`run_method1` und `run_method2` aus Validation.jl auf. Die Adapter verwenden
deine tatsächlichen Verlustfunktionen, Parametergrenzen und Solver; es wird
kein Ersatzoptimierer als validierte Methode ausgegeben.

### Methode 1: Adaptervertrag

Ein NamedTuple `case` enthält:

- `name`: kurzer Dateiname ohne Verzeichnistrenner.
- `initial(rng)`: neuer Parametervektor aus dem übergebenen RNG.
- `optimize(p0)`: Ergebnis `(parameters=..., converged=Bool, residual=...)`.
- `evaluate(p, points)`: Gleichgewicht jeweils neu lösen; Matrix mit Spalten
  Fx [N], Fy [N], Mz [Nm]. Fehler beim Gleichgewichtslösen werfen.
- `target(points)`: Sollwerte mit identischer Spaltenreihenfolge.

Aufruf:

```julia
run_method1(case; seeds=settings["seeds"],
    points=settings["evaluation_points"], directory=output)
```

Die Optimierung selbst verwendet nur optimization_points; evaluate verwendet
zusätzliche Zwischenpunkte. Die CSV enthält komponentenweise MAE, Maximalfehler,
relativen L2-Fehler, Laufzeit und Konvergenz. Bei null Referenznorm ist der relative
Fehler missing. Kräfte und Momente werden nicht ohne Skalierung zusammengerechnet.

### Methode 2: Adaptervertrag

`case` enthält name, edges sowie:

- `initial(rng)`: gleicher Parameterraum und Seeds wie bei den Referenzen.
- `relax(p0)`: `(parameters, beta, converged, residual)` als NamedTuple.
- `admissible(mask)`: prüft mindestens Lageranbindung, Verbindung zum
  Ausgabeknoten und die für den Fall erforderliche kinematische Zulässigkeit.
- `score(parameters, beta)`: `(loss, volume)`, mit identischer dimensionsloser
  Zielfunktion für alle Stufen. Erforderliche Gleichgewichte neu lösen.
- `refine(parameters, mask)`: Methode 1 bei fester diskreter Topologie;
  `(parameters, converged, residual)`.

```julia
run_method2(case; seeds=settings["seeds"], thresholds=settings["thresholds"],
    directory=output, enumerate=true,
    max_edges=settings["max_enumerated_edges"])
```

Pro Seed werden Relaxation, unmittelbare Diskretisierung und Nachoptimierung
separat erfasst. Für sechs Kandidaten werden alle 64 Masken untersucht;
unzulässige und fehlgeschlagene Fälle stehen ebenfalls in der CSV.
Die beste erfolgreich nachoptimierte Referenz ist **kein bewiesenes globales
Optimum**. Der unmittelbare Diskretisierungsschritt übernimmt keinen
Konvergenznachweis der relaxierten Lösung.

Vergleiche loss(discrete)-loss(relaxed) und loss(refined)-loss(discrete)
innerhalb desselben Seeds/Schwellenwerts. Bewerte Balkenanzahl/Volumen gemeinsam
mit Funktionsfehler, Konvergenzquote und Laufzeit. Neue Strukturen müssen das
gleiche Rechenbudget und dieselben Parametergrenzen erhalten.

## Mechanische Tangentenprüfung

`tangent_check(reaction, tangent, q, directions, steps)` vergleicht die
Tangente mit zentralen Differenzen von unabhängig neu gelösten Gleichgewichten.
reaction muss dieselben festgehaltenen/freien DOFs verwenden wie die Reduktion
der Matrix. Verwende mehrere Schrittweiten und konsistent skalierte Richtungen
für Translation und Rotation. Dies überprüft die mechanische Bedeutung von D,
während ForwardDiff lediglich die implementierte Funktion differenziert.

## Ansys-Vergleich

Exportiere Modell und Ansys als numerische CSV mit genau:

```csv
point,Fx,Fy,Mz
```

point ist Verschiebung [m] oder Winkel [rad], für beide Dateien identisch;
Fx/Fy in N, Mz in Nm. Punkte streng aufsteigend, Dezimalpunkt, Komma als Trenner.
Vorzeichen, Ausgabeknoten und Lagerung müssen übereinstimmen. Keine automatische
Interpolation oder Einheitenumrechnung; Modell- und FE-Netzkonvergenz separat dokumentieren.

```sh
julia --project=validation validation/compare_ansys.jl model.csv ansys.csv errors.csv
```

Es werden keine erfundenen Ansys-Daten mitgeliefert. Fehlende FE-Daten gelten
nicht als bestandene Validierung.
