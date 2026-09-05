# mock_data

A small hand-made network for checking the loader and queries without waiting on a real one,
so a failure here is a bug in the loading rather than the chemistry.

- `molecules.csv`, `reactions.csv` — the mock network
- `simple_graph_import.csv` — the same flattened for a direct import
- `example_formatted_data.json` — the shape the loader expects after formatting
- `exported/` — `nodes.txt` and `rels.txt`, the export format
- `test.py` — loads it and runs the checks
