# formose

The formose reaction, the network behind Cruz-Simbron et al., ChemRxiv 2024.

| | |
|---|---|
| `formose_reaction.py` | The run itself. `mod -f formose_reaction.py` |
| `formose_output.txt` | Species with the generation each first appeared in, format in [`main/README.md`](../README.md) |
| `formose_6_may2021.dg` | The generated network dumped by MØD, loadable with `DG.load` without regenerating |
| `formose_rule_count_may2021 (1).tsv` | How many times each rule fired, per generation |
| `OmranDeckerFormoseTestSet.sdf` | Reference structures the generated species were checked against |
| `Formose_Autocatalytic_Cycles_April_2023.zip` | The autocatalytic cycles found in this network, cited by the manuscripts |
| `Neo4j_Imports/` | Exported network, `nodes/` and `rels/`, for loading into Neo4j or Gephi |
