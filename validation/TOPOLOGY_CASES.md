# Referenzfälle für Topologieentscheidungen

`TopologyCases.jl` enthält vier kleine Grundstrukturen, die aus den Optimierungsbeispielen in `mcircle/StructureSynthesis` abgeleitet sind. Jede Struktur besteht ausschließlich aus `Clamp`- und `Branch`-Knoten sowie geraden `Beam`-Elementen. Die Einheiten sind mm, N und Nmm.

| Fall | Knoten | Kandidaten | Referenz | Randbedingung | Geprüfte Entscheidung |
|---|---:|---:|---:|---|---|
| `three_beam_characteristic` | 4 | 6 | 3 | Verschiebung am Clamp, Reaktionskraft als Ziel | Sternstruktur gegen direkte Umgehungskanten |
| `cross_axis_pivot` | 5 | 10 | 4 | Rotation mit zugehöriger Kreisbewegung | gekreuzte Nachgiebigkeiten gegen diagonale Ersatzpfade |
| `snap_through` | 6 | 15 | 5 | horizontale Verschiebung in beide Richtungen | bistabiler Lastpfad gegen diagonale Abkürzung |
| `force_path` | 4 | 6 | 5 | positive und negative Kraft am Branch | Kraftübertragung und Vorzeichen der Ausgangsbewegung |

Die Referenzmaske bezeichnet die aus dem jeweiligen Beispiel übernommene Topologie. Sie ist eine Klassifikationsreferenz für die Kantenentscheidung, kein Beweis, dass diese Topologie unter einer neu gewählten Zielfunktion global optimal ist. Die mechanischen Zielwerte müssen vor der Dissertationsstudie mit dem endgültigen Parametersatz und den gewünschten Kennlinien festgelegt werden.

```julia
include("validation/TopologyCases.jl")
using .TopologyCases

cases = topology_cases(Float32)
case = first(cases)
adj_ground = case.adjacency
adj_reference = topology_adjacency(case)
beams, nodes = prepare_load_case(case, first(case.load_cases))

mask = copy(case.reference)
mask[1] = false
quality = topology_quality(case, mask)
```

`topology_quality` liefert True/False Positives, Precision, Recall, F1, Hamming-Abstand und die graphbasierte Zulässigkeit. `admissible` verlangt eine Anbindung des Ausgabeknotens, die fallspezifischen Mindestgrade der Funktionsknoten und einen Pfad jedes `Branch` zu mindestens einem `Clamp`.

Für die eigentliche Güteauswertung sollten pro Seed und Schwellenwert mindestens folgende Größen gemeinsam berichtet werden:

1. F1 und Hamming-Abstand gegenüber der Referenzmaske,
2. mechanischer Zielfunktionswert nach erneuter Gleichgewichtslösung,
3. Zielfunktionsanstieg durch Diskretisierung,
4. Verbesserung durch Nachoptimierung bei fester Topologie,
5. Volumen beziehungsweise Balkenanzahl und Anteil unzulässiger Topologien.

Eine hohe Übereinstimmung mit der Referenz ist nur dann überzeugend, wenn auch die mechanische Zielfunktion erreicht wird. Umgekehrt kann eine abweichende Topologie gleichwertig oder besser sein; sie sollte deshalb nicht allein wegen einer niedrigeren F1 als falsch verworfen werden.

Die vollständige Enumeration ist für die Fälle mit 6 beziehungsweise 10 Kanten
direkt möglich. Der Snap-through-Fall besitzt 15 Kandidaten; dort sollte
`run_method2(...; enumerate=false)` verwendet oder eine bewusst reduzierte
Kandidatenmenge separat begründet werden.
