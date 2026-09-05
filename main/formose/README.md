# formose

The formose reaction, the network behind Cruz-Simbron et al., ChemRxiv 2024.
Run with `mod -f formose_reaction.py`.

- `formose_reaction.py` — the run
- `formose_output.txt` — species with the generation each appeared in
- `formose_6_may2021.dg` — the network dumped by MØD, loadable with `DG.load`
- `formose_rule_count_may2021 (1).tsv` — how often each rule fired, per generation
- `OmranDeckerFormoseTestSet.sdf` — reference structures checked against
- `Formose_Autocatalytic_Cycles_April_2023.zip` — the cycles found, cited by the manuscripts
- `Neo4j_Imports/` — the network exported for Neo4j or Gephi
