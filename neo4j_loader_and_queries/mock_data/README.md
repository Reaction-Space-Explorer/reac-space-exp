# mock_data

A small hand-made network for checking the loader and the queries without waiting on a real
one. Molecules, reactions and the expected import format, so a failure here is a bug in the
loading rather than in the chemistry.

| | |
|---|---|
| `molecules.csv`, `reactions.csv` | The mock network |
| `simple_graph_import.csv` | The same flattened to one file for a direct import |
| `example_formatted_data.json` | The shape the loader expects after formatting |
| `exported/` | `nodes.txt` and `rels.txt`, the export of this network in the format described in [`main/README.md`](../../main/README.md) |
| `test.py` | Loads the mock data and runs the checks |
