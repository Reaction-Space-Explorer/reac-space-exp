# Neo4j_Imports

The network exported for Neo4j or Gephi, written by `mod_to_neo4j_exporter.py` and split by
generation; the number in each filename is the generation.

- `nodes/nodes_N.txt` — the species: id, SMILES, exact mass, label
- `rels/rels_N.txt` — the reactions: reaction id, SMILES, `-1` reactant or `1` product, rule

`rels_N.txt` is cumulative, so `rels_3.txt` already contains generations 1 and 2. Full
format in [`main/README.md`](../../README.md).
