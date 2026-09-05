# data

Inputs and reference tables used by the generation and filtering steps.

| | |
|---|---|
| `BadRingsList.txt`, `BadAromaticsList.txt` | Forbidden substructures. `main/main.py` reads these to filter species after generation |
| `All_Rules.pdf` | Visual list of the reaction rules in [`rules/`](../rules). May lag the code |
| `Filtered Structures List.pdf` | The structures removed by filtering |
| `DeltaGValuesFinal.xlsx` | Thermochemical data used for the ΔG annotations |
| `formose_filtered_network_deltaG.csv` | Formose network with per-reaction ΔG. Tab-separated in spite of the extension, which is the name the manuscripts cite |
| `*.sdf` | Structure sets: formose test set, HCN, Maillard, alkaline hydrolysis products |
| `*_ms.csv`, `39F.csv`, `2017September21*.csv` | Mass spectra used for the comparisons against generated species |
| `fig7_panel_a.csv`, `fig7_panel_b.csv` | Data behind figure 7 of Arya et al. 2022 |
| `glucose.mol` | Starting structure for the glucose runs |

`All_Rules.pdf`, `Filtered Structures List.pdf` and `formose_filtered_network_deltaG.csv` are
cited by the papers; see [`.github/cited_paths.txt`](../.github/cited_paths.txt).
