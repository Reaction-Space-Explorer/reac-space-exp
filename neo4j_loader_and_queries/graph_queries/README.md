# graph_queries

The Cypher run against the network once it is loaded into Neo4j. This is the pattern match
that finds autocatalytic cycles, and it looks for four structures at once, as the query's own
comments number them:

0. the **ring**,
1. a molecule in the ring splitting to form the starting molecule again, called the
   autocatalytic cord,
2. the **consumers**, drawing species off the ring,
3. the **feeders**, supplying it.

| | |
|---|---|
| `_FINAL_QUERY.TXT` | The query as used |
| `_FINAL_QUERY_PARAMETERIZED.txt` | The same with the ring size left as `{{MIN_RING_SIZE}}` and `{{MAX_RING_SIZE}}`, so one file serves every size rather than being edited each time |
| `example_query.txt` | A worked version with the ring size fixed at 3 to 5, to run against a small network first |
| `_query_visualization_diagram.jpg`, `_query_visualization_diagram2.jpg` | The matched pattern drawn out |
| `node_degree_distribution/` | Degree distribution of the pivot nodes, as a spreadsheet |
| `archive/` | Earlier versions, including a single-label form of the network and a query for the simple loop of Peretó |
