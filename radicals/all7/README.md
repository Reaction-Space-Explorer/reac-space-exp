# all7

A later radical run than [`../all6.py`](../all6.py), adding two things: chains that must not
be formed are filtered out, and the network is exported for Neo4j.

| | |
|---|---|
| `all7.py` | The driver, run with `mod -f all7.py` |
| `forbidden.py` | The chains rejected during generation |
| `rules/` | The 19 rule files this run loads, separated out from the driver |
| `Etapa0.txt` | Starting molecules |
| `data/` | Written by the run: mass density plots and the mass tables behind them, one per generation |
| `Neo4j_Imports/` | Exported network, `nodes/` and `rels/`, in the format described in [`main/README.md`](../../main/README.md) |
