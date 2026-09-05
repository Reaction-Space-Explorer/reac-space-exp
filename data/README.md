# data

Inputs and reference tables for generation and filtering.

- `BadRingsList.txt`, `BadAromaticsList.txt` — forbidden substructures, read by `main/main.py` to filter species after generation
- `All_Rules.pdf` — visual list of the rules in [`rules/`](../rules); may lag the code
- `Filtered Structures List.pdf` — the structures filtering removed
- `DeltaGValuesFinal.xlsx` — thermochemical data behind the ΔG annotations
- `formose_filtered_network_deltaG.csv` — formose network with per-reaction ΔG; tab-separated in spite of the extension, which is the name the manuscripts cite
- `*.sdf` — structure sets: formose test set, HCN, Maillard, alkaline hydrolysis products
- `*_ms.csv`, `39F.csv`, `2017September21*.csv` — mass spectra compared against generated species
- `fig7_panel_a.csv`, `fig7_panel_b.csv` — data behind figure 7 of Arya et al. 2022
- `glucose.mol` — starting structure for the glucose runs
