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
`[0,1]` begrenzt. Eine Gauß-Strafe mit Zentrum 0,5 und der in
`gaussian_sigma` festgelegten Breite drängt Zwischenwerte zu 0 oder 1.

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
binäre Adjazenzmatrix. Methode 2 optimiert kontinuierliche Kantenwerte und
diskretisiert sie anschließend mit dem konfigurierten Schwellwert. Beide Wege
verwenden `Optimisers.Adam`; Iterationszahl, Lernrate, Gleichgewichtsgewicht,
Kraftskala und Schwellwert stehen in `config.toml`. Das Volumen wird weder
berechnet noch bewertet.

Im Ergebnisordner liegen `topology_catalog.csv` sowie je Kennlinie
`*_target.csv`, `*_method1_runs.csv`, `*_topology_summary.csv`,
`*_method2_runs.csv` und `*_method2_comparison.csv`. Bei einem nicht
konvergierten Referenzlauf bleiben die Vergleichswerte leer, statt den gesamten
Versuch abzubrechen.

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
