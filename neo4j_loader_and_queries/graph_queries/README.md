# graph_queries

The Cypher that finds autocatalytic cycles once the network is in Neo4j. It matches four
structures at once, as the query's own comments number them: the ring, a molecule in the ring
splitting to reform the starting molecule (the autocatalytic cord), the consumers, and the
feeders.

- `_FINAL_QUERY.TXT` — the query as used
- `_FINAL_QUERY_PARAMETERIZED.txt` — the same with the ring size left as `{{MIN_RING_SIZE}}`/`{{MAX_RING_SIZE}}`
- `example_query.txt` — a worked version with the ring size fixed at 3 to 5
- `_query_visualization_diagram*.jpg` — the matched pattern drawn out
- `node_degree_distribution/` — degree distribution of the pivot nodes
- `archive/` — earlier versions
