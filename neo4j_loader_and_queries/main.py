import pandas as pd
import numpy as np
import os
from py2neo import Graph, Node, Relationship, NodeMatcher, RelationshipMatcher
import json
import datetime
import matplotlib.pyplot as plt
# from ggplot import *
from shutil import copytree
import math

url = "bolt://neo4j:0000@localhost:7687"
graph = Graph(url)
matcher = NodeMatcher(graph)
rel_matcher = RelationshipMatcher(graph)


def get_timestamp():
    return str(datetime.datetime.now()).replace(":","-").replace(" ","_").replace(".","-")

def rxn_query_str(reactant, product, rxn_id):
    """
    Generate cypher MERGE query for reactant and product node.
    """
    return "MERGE (r: Molecule {smiles_str:\""+ reactant +"\"}) MERGE (p: Molecule {smiles_str:\""+ product +"\"}) MERGE (r)-[:FORMS {rxn_id: "+ str(rxn_id) +"}]->(p)"


def create_reaction_if_not_exists(from_smiles, to_smiles):
    from_molecule = matcher.match("Molecule", smiles_str = from_smiles).first()
    to_molecule = matcher.match("Molecule", smiles_str = to_smiles).first()
    if len(list(graph.match(nodes=(from_molecule, to_molecule), r_type="FORMS"))) <= 0:
        # relationship does not exist
        tx = graph.begin()
        new_r = Relationship(from_molecule, "FORMS", to_molecule)
        tx.create(new_r)
        tx.commit()


def create_relationship_if_not_exists(rxn_id, from_smiles, to_smiles, rule, generation_formed):
    from_molecule = matcher.match("Molecule", smiles_str = from_smiles).first()
    to_molecule = matcher.match("Molecule", smiles_str = to_smiles).first()
    match_pattern = rel_matcher.match(nodes=(from_molecule, to_molecule),
                                r_type="FORMS",
                                # properties = {""}
                                # Remove edge properties so that less edges are formed to reduce
                                # the computational complexity of the queries...
                                # all possible reactions between two molecules can
                                # later be looked up. Properties previously were:
                                properties = {"rule": rule,
                                              "rxn_id": rxn_id,
                                              "generation_formed": generation_formed,
                                              "from_smiles": from_smiles,
                                              "to_smiles": to_smiles
                                              }
                                )
    if len(list(match_pattern)) <= 0:
        # relationship does not exist
        tx = graph.begin()
        # see documentation for weird Relationship function; order of args go:
        # from node, relationship, to node, and then kwargs for relationship properties
        # https://py2neo.org/v4/data.html#py2neo.data.Relationship
        new_r = Relationship(from_molecule, "FORMS", to_molecule, rule=rule, rxn_id=rxn_id, generation_formed = generation_formed)
        tx.create(new_r)
        tx.commit()
        
        # debug (share data for others to import into Neo4j)
        # rels_debug_file = open("mock_data/exported/rels.txt",'a')
        # rels_debug_file.write(f"\n{from_smiles},FORMS,{to_smiles},{rule},{rxn_id},{generation_formed}")
        # rels_debug_file.close()


def create_molecule_if_not_exists(smiles_str, exact_mass, generation_formed):
    """
    Create molecule in DB if not exists.
    """
    molecule = matcher.match("Molecule", smiles_str = smiles_str).first()
    if molecule is None:
        # molecule does not exist, create node with generation information
        tx = graph.begin()
        new_m = Node("Molecule",
                     smiles_str = smiles_str,
                     exact_mass = round(float(exact_mass),3),
                     generation_formed = generation_formed)
        tx.create(new_m)
        tx.commit()
        
        # debug (share data for others to import into Neo4j)
        # nodes_debug_file = open("mock_data/exported/nodes.txt",'a')
        # nodes_debug_file.write(f"\nMolecule,{smiles_str},{exact_mass},{generation_formed}")
        # nodes_debug_file.close()
    else:
        # molecule exists, do nothing
        pass


def import_molecules():
    """
    Takes all output .txt files from data folder and parses the text files for
    any new molecules each generation.
    """
    txt_file_names = os.listdir(os.path.join(os.getcwd(), "data"))
    for file_name in txt_file_names:
        generation = int(file_name.split(".")[0][-1])
        molecules = open(f"data/molecules/{file_name}").read().split("\n")
        for molecule in molecules:
            create_molecule_if_not_exists(smiles_str = molecule,
                                          generation_formed = generation)


def import_reactions():
    pass



def load_simple_graph():
    """
    Connect to Neo4j and load graph.
    
    Usually, graph import is split up into nodes.csv and relationships.csv,
    but since this is a single label & single relationship graph for now,
    this makes the import simpler.
    """
    
    # first, delete all from db
    db.cypher_query("MATCH (n) DETACH DELETE n")
    
    # prepare import data
    data = pd.read_csv("mock_data/simple_graph_import.csv")
    all_molecules = []
    for col in data.columns:
        if col != 'rxn_id':
            all_molecules.extend(list(data[col].unique()))
    all_molecules = list(set(all_molecules)) # get unique set of all molecules and convert back to list
    # all_molecules = [mol for mol in all_molecules if mol != np.nan]
    # load db with import data
    # first, make sure all molecules exist, and create if not with MERGE
    for mol in all_molecules:
        try:
            db.cypher_query("MERGE (:Molecule {smiles_str:\""+mol+"\"})")
        except:
            pass
        
    # then, merge on relationships
    for _, rxn in data.iterrows():
        results, meta1 = db.cypher_query(rxn_query_str(reactant = rxn['reactant_1'], product = rxn['product'], rxn_id = int(rxn['rxn_id'])))
        results, meta2 = db.cypher_query(rxn_query_str(reactant = rxn['reactant_2'], product = rxn['product'], rxn_id = int(rxn['rxn_id'])))
        # print(meta1, meta2) # print meta data on cypher query for each reaction; suppress if too many import records

def import_mock_data():
    """
    Import simple fake dataset to test queries out on.
    """
    # read in data
    molecules = pd.read_csv("mock_data/molecules.csv")
    reactions = pd.read_csv("mock_data/reactions.csv")

    # create nodes
    for _, row in molecules.iterrows():
        create_molecule_if_not_exists(smiles_str = row['smiles_str'],
                                      generation_formed = 0) # doesn't matter yet

    # create relationships
    for _, row in reactions.iterrows():
        create_reaction_if_not_exists(from_smiles = row['from_node'],
                                      to_smiles = row['to_node'])
        # merge_query = rxn_query_str(reactant=row['from_node'],
        #                             product=row['to_node'],
        #                             rxn_id=row['rxn_id'])
        


def import_data_from_MOD_exports(mod_exports_folder_path, generation_limit):
    """
    Clear the database and import the network depending on the Neo4j_Imports
    folder selected.
    """
    # first, clear db
    print("Clearing database...")
    graph.run("MATCH (n) DETACH DELETE n")
    
    # read in data
    nodes_folder = mod_exports_folder_path + "/nodes"
    rels_folder = mod_exports_folder_path + "/rels"

    # create nodes by iterating through each generation's file
    print("Importing nodes...")
    nodes_all_generation_file_names = os.listdir(nodes_folder)
    # first get all generation numbers and sort them in order so the import
    # doesn't go out of order
    all_gens = []
    for generation_file in nodes_all_generation_file_names:
        generation_num = int(generation_file.split("_")[1].split(".")[0])
        all_gens.append(generation_num)
    all_gens.sort()
    
    # set the generation limit to the max (so no data type issue between None/Integer)
    if generation_limit == None:
        generation_limit = max(all_gens)
    
    # create molecule nodes
    for generation_num in all_gens:
        if generation_num <= generation_limit:
            print(f"\tGeneration number {generation_num}...")
            generation_file = "nodes_" + str(generation_num) + ".txt"
            nodes = open(nodes_folder + "/" + generation_file,'r').read().split('\n')
            for node in nodes:
                if node != "":
                    node_data = node.split(',')
                    # print(node_data[0])
                    create_molecule_if_not_exists(smiles_str = node_data[1],
                                                  exact_mass = node_data[2],
                                                  generation_formed = generation_num)

    # create relationships by iterating through each generation's file
    print("Importing relationships...")
    rels_all_generation_file_names = os.listdir(rels_folder)
    # first get all generation numbers and sort them in order so the import
    # doesn't go out of order
    all_rel_gens = []
    for generation_file in rels_all_generation_file_names:
        generation_num = int(generation_file.split("_")[1].split(".")[0])
        all_rel_gens.append(generation_num)
    all_rel_gens.sort()
    for generation_num in all_rel_gens:
        if generation_num <= generation_limit:
            print(f"\tGeneration number {generation_num}...")
            generation_file = "rels_" + str(generation_num) + ".txt"
            rels = open(rels_folder + "/" + generation_file,'r').read().split('\n')
            for rel in rels:
                if rel != "":
                    rel_data = rel.split(',')
                    create_relationship_if_not_exists(rxn_id = rel_data[0],
                                                      from_smiles = rel_data[1],
                                                      to_smiles = rel_data[2],
                                                      rule = rel_data[3],
                                                      generation_formed = generation_num)


def save_query_results(query_result, file_name, this_out_folder):
    with open('output/' + this_out_folder + f"/{file_name}.json", 'w') as file_data_out:
        json.dump(query_result, file_data_out)
    data_df = pd.read_json('output/' + this_out_folder + f"/{file_name}.json")
    data_df.to_csv('output/' + this_out_folder + f"/{file_name}.csv", index=False)

def run_single_value_query(query, value):
    return graph.run(query).data()[0][value]


def get_tabulated_possible_autocatalytic_cycles(mod_exports_folder_path,
                                                this_out_folder,
                                                ring_size_range = (3,7),
                                                feeder_molecule_generation_range = None,
                                                num_structures_limit = 100
                                                ):
    """
    After the graph has been loaded with data, let's execute a query and export
    the tabulated results.
    
    An input of "None" to any of the params means no limit. By default the ring
    size will be from 3 molecules to 7.
    """
    print("\tPreparing query for cycles...")
    
    # make sure inputs are okay
    print("\t\tChecking input parameters...")
    min_ring_size = ring_size_range[0]
    max_ring_size = ring_size_range[1]
    if min_ring_size < 0 or max_ring_size < 0:
        print("Ring sizes can not be negative.")
        quit()
    if min_ring_size > max_ring_size:
        print("The minimum ring size must not exceed the maximum.")
        quit()
    if min_ring_size <= 2:
        print("The minimum ring size must be above 2.")
        quit()
    
    if feeder_molecule_generation_range != None:
        min_feeder_gen = feeder_molecule_generation_range[0]
        max_feeder_gen = feeder_molecule_generation_range[1]
        if min_feeder_gen < 0 or max_feeder_gen < 0:
            print("The feeder generation can not be negative.")
            quit()
        if min_feeder_gen > max_feeder_gen:
            print("The minimum feeder generation must not exceed the maximum.")
            quit()
    else:
        min_feeder_gen = None
        max_feeder_gen = None
    
    # load query and insert params
    print("\t\tReplacing query parameters in query string...")
    query_txt = open("graph_queries/_FINAL_QUERY_PARAMETERIZED.txt",'r').read()
    query_txt = query_txt.replace("{{MIN_RING_SIZE}}", str(min_ring_size))
    query_txt = query_txt.replace("{{MAX_RING_SIZE}}", str(max_ring_size))
    
    
    if feeder_molecule_generation_range == None:
        query_txt = query_txt.replace("{{COMMENT_OUT_FEEDER_GEN_LOGIC}}", "//")
    else:
        query_txt = query_txt.replace("{{COMMENT_OUT_FEEDER_GEN_LOGIC}}", "")
        query_txt = query_txt.replace("{{MIN_FEEDER_GENERATION}}", str(min_feeder_gen))
        query_txt = query_txt.replace("{{MAX_FEEDER_GENERATION}}", str(max_feeder_gen))
    
    query_txt = query_txt.replace("{{NUM_STRUCTURES_LIMIT}}", str(num_structures_limit))
    
    # Get the max ID of all molecules to get a random molecule to start with.
    # Query several times in small chunks to stochastically estimate the behavior
    # of the graph without having to traverse the entire thing for this query.
    # max_node_id = run_single_value_query("MATCH (n) RETURN max(ID(n)) AS max_node_id","max_node_id")
    # WHERE ID(beginMol) = round(rand() * {{MAX_NODE_ID}})
    # print("\t\t\t" + query_txt)
    
    # Execute query in Neo4j. If out of memory error occurs, need to change DB settings:
    # I used heap initial size set to 20G, heap max size set to 20G, and page cache size set to 20G,
    # but these settings would depend on your hardware limitations.
    # See Neo4j Aura for cloud hosting: https://neo4j.com/aura/
    print("\t\tExecuting query and collecting results (this may take awhile)...")
    print(f"\t\t\tTime start: {get_timestamp()}")
    query_result = graph.run(query_txt).data()
    print(f"\t\t\tTime finish: {get_timestamp()}")
    # print("\t\tQuery results:")
    # print(query_result[0])
    print("\t\tSaving query results and meta info...")
    
    # save data as JSON and CSV (JSON for easy IO, CSV for human readability)
    save_query_results(query_result = query_result,
                       file_name = "autocat_query_results",
                       this_out_folder = this_out_folder)
    
    # save meta info as well in out folder
    with open("output/" + this_out_folder + "/autocat_query.txt", 'w') as file_query_out:
        file_query_out.write(query_txt)
    query_params = pd.DataFrame( {"parameter": ["min_ring_size","max_ring_size","min_feeder_gen","max_feeder_gen","num_structures_limit"],
                                  "value": [min_ring_size, max_ring_size, min_feeder_gen, max_feeder_gen, num_structures_limit] } )
    query_params.to_csv("output/" + this_out_folder + "/autocat_query_parameters.csv", index=False)
    return this_out_folder



def analyze_possible_autocatalytic_cycles(mod_exports_folder_path, query_results_folder):
    """
    Now that we have the tabulated results of the graph queries, let's do some
    analysis on what's going on.
    
    1. Ring size frequency distribution
    2. Total mass per cycle per feeder molecule's generation (calculate total
        using only the molecules in the ring, and use the feeder molecule's
        generation as the ring's generation).
        Note: make sure to remove duplicates when getting sum of mass in ringPathNodes
        because the beginMol is counted twice (it is the start and end node in the path).
    3. Count of cycles by feeder generation
    """
    
    
    print("Generating some plots on cycle size distribution / stats by generation...")
    # 1.
    query_data = pd.read_json("output/" + query_results_folder + "/autocat_query_results.json")
    # print(query_data.describe())
    # print(query_data.head())
    
    # cycle distribution (y axis is frequency, x axis is ring size)
    fig, ax = plt.subplots()
    query_data['countMolsInRing'].value_counts().plot(ax = ax,
                                                      kind='bar',
                                                      title = "Ring Size Frequency Distribution")
    ax.set_xlabel("Ring Size (# of Molecules)")
    ax.set_ylabel("Count of Cycles")
    plt.savefig("output/" + query_results_folder + "/ring_size_distribution.png")
    # plt.show()
    
    
    # 2.
    # Total mass of cycle per generation. Not really needed.
    
    
    # 3.
    # count of cycles by feeder generation
    fig, ax = plt.subplots()
    gen_formed_arr = []
    feederMolData = list(query_data['feederMol'])
    for feederMol in feederMolData:
        gen_formed_arr.append(feederMol['generation_formed'])
    # get unique list of feeder generations and sum by generation
    gen_formed_arr = np.array(gen_formed_arr)
    feeder_gen_counts = np.unique(gen_formed_arr, return_counts=True)
    feeder_gen_counts = np.transpose(feeder_gen_counts)
    cycles_by_gen_df = pd.DataFrame(feeder_gen_counts, columns=['feeder_gen',
                                                                'cycle_count'])
    cycles_by_gen_df.plot(ax=ax,
                          x = "feeder_gen",
                          y = "cycle_count",
                          kind = "bar",
                          legend = False,
                          title = "Count of Cycles by Feeder Generation")
    ax.set_xlabel("Cycle Generation (Generation Formed of Feeder Molecule)")
    ax.set_ylabel("Count of Cycles")
    plt.savefig("output/" + query_results_folder + "/count_cycles_by_feeder_generation.png")
    
    print("\tAutocatalysis pattern matching done.")




def plot_hist(file_name, statistic_col_name, title, x_label, y_label):
    # fig, ax = plt.subplots()
    # df = pd.read_csv(f"output/{query_results_folder}/{file_name}.csv")
    # num_bins = int(math.sqrt(df.shape[0])) # estimate the number of bins by taking the square root of the number of rows in the dataset
    # df.plot.hist(bins=num_bins, ax=ax)
    # ax.set_xlabel(x_label)
    # ax.set_ylabel(y_label)
    # plt.savefig(f"output/{query_results_folder}/{file_name}.png")
    fig, ax = plt.subplots()
    df = pd.read_csv(f"output/{query_results_folder}/{file_name}.csv")
    num_bins = int(math.sqrt(df.shape[0])) # estimate the number of bins by taking the square root
    df = pd.pivot_table(df,
                        values="smiles_str",
                        index=[statistic_col_name],
                        columns=["generation_formed"],
                        aggfunc=lambda x: math.sqrt(len(x.unique()))) # the square of the count of unique smiles_str
    df.plot.hist(ax=ax,
                 bins = num_bins,
                 title=title,
                 figsize = (15,15))
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    plt.savefig(f"output/{query_results_folder}/{file_name}_histogram.png")


def plot_scatter(file_name, statistic_col_name, title, x_label, y_label):
    fig, ax = plt.subplots()
    df = pd.read_csv(f"output/{query_results_folder}/{file_name}.csv")
    df = df.head(100) # cut off by top 100 most interesting 
    # df.plot.bar(ax = ax,
    #             x = "smiles_str",
    #             y = statistic_col_name,
    #             # color = "generation_formed",
    #             legend = True,
    #             title = title,
    #             figsize = (14,14))
    # ggplot(aes(x = "smiles_str",
    #            y = statistic_col_name,
    #            color = "generation_formed"),
    #        data = df) + geom_point()
    # ax.legend(['generation_formed'])
    groups = df.groupby("generation_formed")
    for name, group in groups:
        plt.plot(group['smiles_str'],
                 group[statistic_col_name],
                 marker = "o",
                 linestyle = "",
                 label = name)
    plt.legend(loc="best", title="Generation Formed")
    plt.xticks(rotation=90)
    fig.set_figheight(15)
    fig.set_figwidth(15)
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    plt.savefig(f"output/{query_results_folder}/{file_name}_scatter.png")


def network_statistics(query_results_folder):
    """
    Get some statistics on the network.
    0. Number of nodes and edges in the graph, as well as various network-level
        statistics: 1. Eigenvector centrality, 2. Betweenness Centrality,
        3. Random-walk betweenness, 4. Clique enumeration,
        5. k-plex enumeration, 6. k-core enumeration,
        7. k-component enumeration, 8. neighbor redundancy
    1. Node degree distribution: sqrt of node degree frequency by degree
        value colored by generation_formed, one plot for incoming, outgoing,
        and incoming and outgoing edges
    2. Avg number of edges per node per generation
    """
    
    print("Doing some network statistics...")
    # 0.
    # get total number of nodes and edges
    total_count_nodes = run_single_value_query("MATCH (n) RETURN COUNT(n) AS count_nodes", 'count_nodes')
    total_count_rels = run_single_value_query("MATCH (n)-[r]->() RETURN COUNT(r) AS count_rels", 'count_rels')
    
    # 0.1 eigenvector_centrality
    # do by generation, molecule, order by score first
    eigenvector_centrality = graph.run("""
                                       CALL algo.eigenvector.stream('Molecule', 'FORMS', {})
                                       YIELD nodeId, score
                                       RETURN algo.asNode(nodeId).smiles_str AS smiles_str, algo.asNode(nodeId).generation_formed AS generation_formed, score
                                       ORDER BY score DESC """).data()
    save_query_results(eigenvector_centrality, "eigenvector_centrality", query_results_folder)
    plot_hist(file_name = "eigenvector_centrality",
              statistic_col_name = "score",
              title = "Histogram of Eigenvector Centrality",
              x_label = "Eigenvector Centrality Score Bin",
              y_label = "Count of Molecules")
    plot_scatter(file_name = "eigenvector_centrality",
                 statistic_col_name = "score",
                 title = "Eigenvector Centrality - Top 100 Connected Molecules",
                 x_label = "Molecule Smiles Format",
                 y_label = "Eigenvector Centrality Score")
    avg_eigenvector_centrality = run_single_value_query("""
                                                        CALL algo.eigenvector.stream('Molecule', 'FORMS', {})
                                                        YIELD nodeId, score
                                                        RETURN avg(score) AS avg_score
                                                        """,
                                                        "avg_score")
    
    
    # 0.2 betweenness_centrality
    betweenness_centrality = graph.run("""
                                       CALL algo.betweenness.stream('Molecule','FORMS',{direction:'out'})
                                       YIELD nodeId, centrality
                                       MATCH (molecule:Molecule) WHERE id(molecule) = nodeId
                                       RETURN molecule.smiles_str AS smiles_str, molecule.generation_formed AS generation_formed, centrality
                                       ORDER BY centrality DESC;
                                       """).data()
    save_query_results(betweenness_centrality, "betweenness_centrality", query_results_folder)
    plot_hist(file_name = "betweenness_centrality",
              statistic_col_name = "centrality",
              title = "Histogram of Betweenness Centrality",
              x_label = "Betweenness Centrality Score Bin",
              y_label = "Count of Molecules")
    plot_scatter(file_name = "betweenness_centrality",
                 statistic_col_name = "centrality",
                 title = "Betweenness Centrality - Top 100 Connected Molecules",
                 x_label = "Molecule Smiles Format",
                 y_label = "Betweenness Centrality Score")
    avg_betweenness_centrality = run_single_value_query("""
                                       CALL algo.betweenness.stream('Molecule','FORMS',{direction:'out'})
                                       YIELD nodeId, centrality
                                       MATCH (molecule:Molecule) WHERE id(molecule) = nodeId
                                       RETURN avg(centrality) AS avg_centrality
                                       """, 'avg_centrality')
    
    # 0.3 Random-walk betweenness
    random_walk_betweenness = graph.run(""" CALL algo.betweenness.sampled.stream('Molecule','FORMS', {strategy:'random', probability:1.0, maxDepth:1, direction: "out"})
                                        YIELD nodeId, centrality
                                        MATCH (molecule) WHERE id(molecule) = nodeId
                                        RETURN molecule.smiles_str AS smiles_str, molecule.generation_formed AS generation_formed, centrality AS random_walk_centrality
                                        ORDER BY random_walk_centrality DESC;""").data()
    save_query_results(random_walk_betweenness, "random_walk_betweenness", query_results_folder)
    plot_hist(file_name = "random_walk_betweenness",
              statistic_col_name = "random_walk_centrality",
              title = "Histogram of Random Walk Betweenness Centrality",
              x_label = "Random Walk Betweenness Centrality Score Bin",
              y_label = "Count of Molecules")
    plot_scatter(file_name = "random_walk_betweenness",
                 statistic_col_name = "random_walk_centrality",
                 title = "Random Walk Betweenness Centrality - Top 100 Connected Molecules",
                 x_label = "Molecule Smiles Format",
                 y_label = "Random Walk Betweenness Centrality Score")
    avg_random_walk_betweenness = run_single_value_query("""CALL algo.betweenness.stream('Molecule','FORMS',{direction:'out'})
                                                         YIELD nodeId, centrality
                                                         MATCH (molecule:Molecule) WHERE id(molecule) = nodeId
                                                         RETURN avg(centrality) AS avg_random_walk_betweenness""",
                                                         'avg_random_walk_betweenness')
    
    # 0.4 Clique enumeration
    avg_clique_enumeration = None #run_single_value_query("", 'clique_enumeration')
    
    # 0.5 K-Plex enumeration
    avg_k_plex_enumeration = None #run_single_value_query("", 'k_plex_enumeration')
    
    # 0.6 K-Core enumeration
    avg_k_core_enumeration = None #run_single_value_query("", 'k_core_enumeration')
    
    # 0.7 K-Component enumeration
    avg_k_component_enumeration = None #run_single_value_query("", 'k_component_enumeration')
    
    # 0.8 Neighbor redundancy
    avg_neighbor_redundancy = None #run_single_value_query("", 'neighbor_redundancy')
    
    # save all to graph_info DataFrame
    graph_info = pd.DataFrame({"statistic": ["Total Count Molecules", "Total Count Edges","Average Eigenvector Centrality", "Average Betweenness Centrality", "Average Random-walk Betweenness", 'Clique enumeration','k-plex enumation','k-core enumeration','k-component enumeration','Neighbor redundancy'],
                               "value": [total_count_nodes, total_count_rels, avg_eigenvector_centrality, avg_betweenness_centrality, avg_random_walk_betweenness, avg_clique_enumeration, avg_k_plex_enumeration, avg_k_core_enumeration, avg_k_component_enumeration, avg_neighbor_redundancy]})
    graph_info.to_csv(f"output/{query_results_folder}/_network_info.csv", index=False)
    
    
    # 1.
    # first do the query and save the results
    node_deg_query = """
    MATCH (n:Molecule)
    RETURN n.smiles_str AS smiles_str, n.generation_formed AS generation_formed, size((n)--()) AS count_relationships
    """
    node_deg_query_results = graph.run(node_deg_query).data()
    node_deg_file = "node_distribution_results"
    save_query_results(query_result = node_deg_query_results,
                       file_name = node_deg_file,
                       this_out_folder = query_results_folder)
    
    # now read in the results, transform, and plot
    # also can represent this as a histogram?
    fig, ax = plt.subplots()
    node_deg_df = pd.read_csv(f"output/{query_results_folder}/{node_deg_file}.csv")
    # node_deg_df['count_relationships'].value_counts().plot(ax=ax,
    #                                                        kind='bar',
    #                                                        title="Node Degree Distribution by Generation Formed")
    node_deg_pivot = pd.pivot_table(node_deg_df,
                                    values="smiles_str",
                                    index=["count_relationships"],
                                    columns=["generation_formed"],
                                    aggfunc=lambda x: math.sqrt(len(x.unique()))) # the square of the count of unique smiles_str
    node_deg_pivot.plot(ax=ax,
                        kind="bar",
                        title="Square of Molecule Degree by Generation Formed",
                        figsize = (8,5))
    ax.set_xlabel("Molecule Degree (count of incoming and outgoing edges)")
    ax.set_ylabel("sqrt(Count of Molecules)")
    plt.savefig(f"output/{query_results_folder}/{node_deg_file}.png")
    
    # 2.
    # get average number of edges by node and generation
    fig, ax = plt.subplots()
    node_deg_avg = node_deg_df.groupby(by=['generation_formed']).mean().reset_index()
    # print(node_deg_avg)
    node_deg_avg.plot(ax=ax,
                      x = "generation_formed",
                      y = "count_relationships",
                      kind="scatter",
                      title="Average Molecule Degree by Generation Formed",
                      figsize = (8,5),
                      legend = False)
    ax.set_xlabel("Generation Formed")
    ax.set_ylabel("Average Node Degree")
    plt.savefig(f"output/{query_results_folder}/{node_deg_file}_avg.png")
    
    # incoming relationships by molecule
    incoming_rels_count_file = "incoming_rels_count"
    incoming_rels_count = graph.run("""
                                    MATCH (n)<-[r:FORMS]-()
                                    RETURN n.smiles_str AS smiles_str,
                                    n.generation_formed AS generation_formed,
                                    n.exact_mass AS exact_mass,
                                    count(r) AS count_incoming
                                    ORDER BY count_incoming DESC
                                    """).data()
    save_query_results(query_result = incoming_rels_count,
                       file_name = incoming_rels_count_file,
                       this_out_folder = query_results_folder)
    fig, ax = plt.subplots()
    node_deg_df = pd.read_csv(f"output/{query_results_folder}/{incoming_rels_count_file}.csv")
    # node_deg_df['count_relationships'].value_counts().plot(ax=ax,
    #                                                        kind='bar',
    #                                                        title="Node Degree Distribution by Generation Formed")
    node_deg_pivot = pd.pivot_table(node_deg_df,
                                    values="smiles_str",
                                    index=["count_incoming"],
                                    columns=["generation_formed"],
                                    aggfunc=lambda x: math.sqrt(len(x.unique()))) # the square of the count of unique smiles_str
    node_deg_pivot.plot(ax=ax,
                        kind="bar",
                        title="Square of Molecule Degree by Generation Formed for Incoming Relationships",
                        figsize = (8,5))
    ax.set_xlabel("Molecule Degree (count of incoming edges)")
    ax.set_ylabel("sqrt(Count of Molecules)")
    plt.savefig(f"output/{query_results_folder}/{incoming_rels_count_file}.png")
    
    # outgoing relationships by molecule
    outgoing_rels_count_file = "outgoing_rels_count"
    outgoing_rels_count = graph.run("""
                                    MATCH (n)-[r:FORMS]->()
                                    RETURN n.smiles_str AS smiles_str,
                                    n.generation_formed AS generation_formed,
                                    n.exact_mass AS exact_mass,
                                    count(r) AS count_outgoing
                                    ORDER BY count_outgoing DESC
                                    """).data()
    save_query_results(query_result = outgoing_rels_count,
                       file_name = outgoing_rels_count_file,
                       this_out_folder = query_results_folder)
    fig, ax = plt.subplots()
    node_deg_df = pd.read_csv(f"output/{query_results_folder}/{outgoing_rels_count_file}.csv")
    # node_deg_df['count_relationships'].value_counts().plot(ax=ax,
    #                                                        kind='bar',
    #                                                        title="Node Degree Distribution by Generation Formed")
    node_deg_pivot = pd.pivot_table(node_deg_df,
                                    values="smiles_str",
                                    index=["count_outgoing"],
                                    columns=["generation_formed"],
                                    aggfunc=lambda x: math.sqrt(len(x.unique()))) # the square of the count of unique smiles_str
    node_deg_pivot.plot(ax=ax,
                        kind="bar",
                        title="Square of Molecule Degree by Generation Formed for Outgoing Relationships",
                        figsize = (8,5))
    ax.set_xlabel("Molecule Degree (count of outgoing edges)")
    ax.set_ylabel("sqrt(Count of Molecules)")
    plt.savefig(f"output/{query_results_folder}/{outgoing_rels_count_file}.png")
    
    
    
    
    print("\tNetwork statistics done.")



if __name__ == "__main__":
    # choose a path for the Neo4j_Imports folder to import the data from MOD into Neo4j
    mod_exports_folder_path = "../main/Neo4j_Imports"
    # mod_exports_folder_path = "../radicals/all7/Neo4j_Imports"
    # import_data_from_MOD_exports(mod_exports_folder_path = mod_exports_folder_path,
    #                              generation_limit = 2) # Set to None or Integer. The generation limit at which to import
    
    # create a timestamped output folder to store everything for this run
    query_results_folder = get_timestamp()
    os.mkdir("output/" + query_results_folder)
    # query_results_folder = "2020-08-01_18-39-23-986232" # manually override folder name for debugging
    
    # Optional: save the state of the Neo4j_Imports folder at the time this was run
    copytree(mod_exports_folder_path, 'output/' + query_results_folder + "/Neo4j_Imports")
    
    # do pattern match query on possible autocatalytic cycles
    # get_tabulated_possible_autocatalytic_cycles(mod_exports_folder_path = mod_exports_folder_path,
    #                                             this_out_folder = query_results_folder,
    #                                             ring_size_range = (3, 5),
    #                                             feeder_molecule_generation_range = None,
    #                                             num_structures_limit = 100) #75000
    # analyze_possible_autocatalytic_cycles(mod_exports_folder_path = mod_exports_folder_path,
    #                                       query_results_folder = query_results_folder)
    
    # do network statistics and get plots
    network_statistics(query_results_folder = query_results_folder)














