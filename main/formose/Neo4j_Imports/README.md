# Neo4j_Imports

The generated network exported for loading into Neo4j, or into Gephi for drawing. Written by
`mod_to_neo4j_exporter.py`, and split by generation: the number in each filename is the
generation it covers.

| | |
|---|---|
| `nodes/nodes_N.txt` | The species. Tab separated: id, SMILES, exact mass, and the label `Molecule` |
| `rels/rels_N.txt` | The reactions. Tab separated: reaction id, molecule SMILES, `-1` for a reactant or `1` for a product, and the rule that fired |

`rels_N.txt` is cumulative, so `rels_3.txt` already contains generations 1 and 2. Getting one
generation alone means subtracting the previous file. The full format, including why reaction
ids carry a subscript, is in [`main/README.md`](../../README.md).
