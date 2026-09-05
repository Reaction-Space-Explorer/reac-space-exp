# all7

A later radical run than [`../all6.py`](../all6.py), filtering out chains that must not form
and exporting for Neo4j. Run with `mod -f all7.py`.

- `all7.py` — the driver
- `forbidden.py` — the chains rejected during generation
- `rules/` — the 19 rules this run loads
- `Etapa0.txt` — starting molecules
- `data/` — mass density plots and the mass tables behind them, one per generation
- `Neo4j_Imports/` — the network exported for Neo4j
