This output folder is where the data from each run is collected into folder named by the timestamp it was created.

Below is an explanation of what each file in these folders means:

1. Neo4j_Imports Folder:
    This folder is copied from the source Neo4j_Imports folder and recorded here so that if you re-generate the reaction network with different parameters, you still have the network data from MOD from which this Neo4j network was generated.

2. query.txt:
    This file is the exact string passed to Neo4j to execute. The query parameters were passed programmatically, by the function call get_tabulated_possible_autocatalytic_cycles(). See also the file ../graph_queries/_FINAL_QUERY_PARAMETERIZED.txt. Query parameters are found and replaced using {{VARIABLE}} notation.

3. query_parameters.csv:
    This is a simple recording of the variables passed through get_tabulated_possible_autocatalytic_cycles() at the time the script was ran.

4. query_results.csv:
    This is the raw data output of the query, where each row represents a pattern matched by the query (where our target pattern is possible autocatalytic cycles). Has the same exact data as query_results.json, only it is more "human readable". To understand what the structure names in the data mean, refer to ../graph_queries/_query_visualization_diagram.jpg.

5. query_results.json:
    This is the raw data output of the query, where each index in the array represents a pattern matched by the query (where our target pattern is possible autocatalytic cycles). Has the same exact data as query_results.csv, only this format is easier to convert between dictionary/array objects. To understand what the structure names in the data mean, refer to ../graph_queries/_query_visualization_diagram.jpg.

6. ring_size_distribution.png:
    This plot is a simple plot of the number of ring patterns found by the number of molecules in the ring. Note that this is highly unreliable when matching a small number of patterns.



