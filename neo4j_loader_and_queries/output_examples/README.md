# output_examples

One complete run of the loader and queries, kept so the output format can be seen without
running anything. The folder is named for when it ran.

Inside `2020-08-21_19-20-09-619330/`:

| | |
|---|---|
| `autocat_query.txt`, `autocat_query_parameters.csv` | The query run, and the parameters it was given |
| `autocat_query_results.csv`, `.json` | The cycles it matched |
| `betweenness_centrality.csv`, `.json` | Centrality over the network |
| `_network_info.csv` | Size and shape of the network queried |
| `Neo4j_Imports/` | The network itself, in the import format |
