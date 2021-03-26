import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
import os
from py2neo import Graph, Node, Relationship, NodeMatcher, RelationshipMatcher
# from neo4j import GraphDatabase
# import neo4j
import networkx as nx
import json
import datetime
import matplotlib.pyplot as plt
# from ggplot import *
from shutil import copytree
import math
# from graph_tool.all import *
import json
import random


# Choose a path for the Neo4j_Imports folder to import the data from MOD into Neo4j
# formose_MOD_exports_path = "../data/formose/Neo4j_Imports"
formose_MOD_exports_path = "../data/pyruvic_acid/Neo4j_Imports"
glucose_MOD_exports_path = "../data/glucose/Neo4j_Imports"
# exports_folder_paths = [formose_MOD_exports_path, glucose_MOD_exports_path]
EXPORT_PATHS = [glucose_MOD_exports_path]

# Set the following to False if you want to leave order of import records in
# each generation file the same; set to True to randomly shuffle the order of
# the records within each file. By shuffling the order, the order at which the
# molecules are imported into Neo4j will be randomized, and thus the start point
# at which the cycles pattern match begins is randomized each time, so we can
# get samples at different starting points in the network since it is too
# computationally intensive to match for all possible patterns in the network.
SHUFFLE_GENERATION_DATA = True
# Repeat the whole import and pattern match routine REPEAT_RUNS amount of times.
# Pair this with SHUFFLE_GENERATION_DATA so that if SHUFFLE_GENERATION_DATA
# is True, sample pattern matches on the graph REPEAT_RUNS amount of times
# starting from random points on the graph from the shuffling, where each
# run matches up to NUM_STRUCTURES_LIMIT of patterns.
REPEAT_RUNS = 10

# Filter out these molecules by smiles string from being imported into Neo4j
# for pattern match / network statistic calculations.
MOLECULE_FILTER = ['O']

# If True, will match for autocatalytic pattern mattches using the pattern match
# query in graph_queries/_FINAL_QUERY_PARAMETERIZED.txt. If not, will skip this
# and just do node degree / rank calculations. (One reason you might want to disable
# pattern match query results is because this is very computationally intensive
# and takes a lot of time; so disable if you are just looking for network statistics.)
PATTERN_MATCHES = True
# Rather than disabling completely if running into performance issues, limit the
# number of patterns that can be matched so that the query stops executing as
# soon as it reaches the pattern limit, and the matches are returned.
NUM_STRUCTURES_LIMIT = 100

# Limit the number of generations that each network can be imported on. If None,
# no limit--will default to the maximum number of generations generated. You may
# want to limit this to ~4 generations or less if performance is an issue; the
# network will grow exponentially, so pattern match queries might take too long
# to produce results.
GENERATION_LIMIT = 4 # None

# If NETWORK_SNAPSHOTS is True, the program gathers data on the network at each generation
# in the reaction netowrk. If False, the program gathers data only on the state of
# the network once all generations have completely finished being loaded (snapshot
# only of the final generation).
NETWORK_SNAPSHOTS = True

# Enable this only if you want to capture network statistics (such as node degree
# plots over generation)
COLLECT_NETWORK_STATISTICS = False

# Set this to True if you want to generate a static image of the network after
# loading. Might run into Out of Memory error. Default leaving this as False
# because we generated a much nicer visualization of the full network using Gephi.
FULL_NETWORK_VISUALIZATION = False


# configure network database Neo4j
url = "bolt://neo4j:0000@localhost:7687"
graph = Graph(url)
matcher = NodeMatcher(graph)
rel_matcher = RelationshipMatcher(graph)


def get_timestamp():
    return str(datetime.datetime.now()).replace(":","-").replace(" ","_").replace(".","-")


def create_molecule_if_not_exists(smiles_str, generation_formed, exact_mass=0):
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
        return new_m
    return molecule

def create_reaction_if_not_exists(id, rule, generation_formed):
    reaction = matcher.match("Reaction", id = id).first()
    if reaction is None:
        tx = graph.begin()
        new_rxn = Node("Reaction",
                       id = id,
                       rule = rule,
                       generation_formed = generation_formed)
        tx.create(new_rxn)
        tx.commit()
        return new_rxn
    return reaction
    

def create_reactant_rel_if_not_exists(smiles_str, rxn_id, generation_formed):
    molecule = matcher.match("Molecule", smiles_str = smiles_str).first()
    reaction = matcher.match("Reaction", id = rxn_id).first()
    match_pattern = rel_matcher.match(nodes=(molecule, reaction),
                                      r_type="REACTANT" #,
                                      # properties = {"generation_formed": generation_formed}
                                      )
    # if pattern does not exist in db
    if len(list(match_pattern)) <= 0:
        tx = graph.begin()
        # see documentation for weird Relationship function; order of args go:
        # from node, relationship, to node, and then kwargs for relationship properties
        # https://py2neo.org/v4/data.html#py2neo.data.Relationship
        new_r = Relationship(molecule, "REACTANT", reaction,
                             generation_formed=generation_formed)
        tx.create(new_r)
        tx.commit()
        return new_r
    return match_pattern

def create_product_rel_if_not_exists(smiles_str, rxn_id, generation_formed):
    molecule = matcher.match("Molecule", smiles_str = smiles_str).first()
    reaction = matcher.match("Reaction", id = rxn_id).first()
    match_pattern = rel_matcher.match(nodes=(reaction, molecule),
                                      r_type="PRODUCT" #,
                                      # properties = {"generation_formed": generation_formed}
                                      )
    # if pattern does not exist in db
    if len(list(match_pattern)) <= 0:
        tx = graph.begin()
        # see documentation for weird Relationship function; order of args go:
        # from node, relationship, to node, and then kwargs for relationship properties
        # https://py2neo.org/v4/data.html#py2neo.data.Relationship
        new_r = Relationship(reaction, "PRODUCT", molecule,
                             generation_formed=generation_formed)
        tx.create(new_r)
        tx.commit()
        return new_r
    return match_pattern


def save_query_results(generation_num, query_result, file_name, this_out_folder):
    with open(f'output/' + this_out_folder + f"/{generation_num}/{file_name}.json", 'w') as file_data_out:
        json.dump(query_result, file_data_out)
    data_df = pd.read_json(f'output/' + this_out_folder + f"/{generation_num}/{file_name}.json")
    data_df.to_csv(f'output/' + this_out_folder + f"/{generation_num}/{file_name}.csv", index=False)


def read_query_results(file_path):
    try:
        df = pd.read_csv(file_path)
    except:
        df = pd.DataFrame()
    return df


def run_single_value_query(query, value):
    return graph.run(query).data()[0][value]


def get_tabulated_possible_autocatalytic_cycles(generation_num,
                                                mod_exports_folder_path,
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
    print("\t\t\tPreparing query for cycles...")
    
    # make sure inputs are okay
    print("\t\t\t\tChecking input parameters...")
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
    print("\t\t\t\tReplacing query parameters in query string...")
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
    print("\t\t\t\tExecuting query and collecting results (this may take awhile)...")
    print(f"\t\t\t\tTime start: {get_timestamp()}")
    query_result = graph.run(query_txt).data()
    print(f"\t\t\t\tTime finish: {get_timestamp()}")
    # print("\t\tQuery results:")
    # print(query_result[0])
    print("\t\t\t\tSaving query results and meta info...")
    
    # save data as JSON and CSV (JSON for easy IO, CSV for human readability)
    save_query_results(generation_num = generation_num,
                       query_result = query_result,
                       file_name = "autocat_query_results",
                       this_out_folder = this_out_folder)
    
    # save meta info as well in out folder
    with open(f"output/" + this_out_folder + f"/{generation_num}/autocat_query.txt", 'w') as file_query_out:
        file_query_out.write(query_txt)
    query_params = pd.DataFrame( {"parameter": ["min_ring_size","max_ring_size","min_feeder_gen","max_feeder_gen","num_structures_limit"],
                                  "value": [min_ring_size, max_ring_size, min_feeder_gen, max_feeder_gen, num_structures_limit] } )
    query_params.to_csv(f"output/" + this_out_folder + f"/{generation_num}/autocat_query_parameters.csv", index=False)
    return this_out_folder



def analyze_possible_autocatalytic_cycles(generation_num, mod_exports_folder_path, query_results_folder):
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
    query_data = pd.read_json(f"output/" + query_results_folder + f"/{generation_num}/autocat_query_results.json")
    if query_data.empty:
        print("No cycles found.")
        return
    # print(query_data.describe())
    # print(query_data.head())
    
    # cycle distribution (y axis is frequency, x axis is ring size)
    fig, ax = plt.subplots()
    # print(query_data.head())
    # query_data['countMolsInRing'] = query_data['countMolsInRing'].astype(int)
    query_data['countMolsInRing'].value_counts().plot(ax = ax,
                                                      kind='bar',
                                                      title = "Ring Size Frequency Distribution")
    ax.set_xlabel("Ring Size (# of Molecules)")
    ax.set_ylabel("Count of Cycles")
    plt.savefig(f"output/" + query_results_folder + f"/{generation_num}/ring_size_distribution.png")
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
    plt.savefig(f"output/" + query_results_folder + f"/{generation_num}/count_cycles_by_feeder_generation.png")
    
    # close all plots so they don't accumulate memory
    print("\tAutocatalysis pattern matching done.")
    plt.close('all')




def plot_hist(query_results_folder, generation_num, file_name, statistic_col_name, title, x_label, y_label):
    # fig, ax = plt.subplots()
    # df = pd.read_csv(f"output/{query_results_folder}/{file_name}.csv")
    # num_bins = int(math.sqrt(df.shape[0])) # estimate the number of bins by taking the square root of the number of rows in the dataset
    # df.plot.hist(bins=num_bins, ax=ax)
    # ax.set_xlabel(x_label)
    # ax.set_ylabel(y_label)
    # plt.savefig(f"output/{query_results_folder}/{file_name}.png")
    fig, ax = plt.subplots()
    df = pd.read_csv(f"output/{query_results_folder}/{generation_num}/{file_name}.csv")
    num_bins = int(math.sqrt(df.shape[0])) # estimate the number of bins by taking the square root
    df = pd.pivot_table(df,
                        values="smiles_str",
                        index=[statistic_col_name],
                        columns=["generation_formed"],
                        aggfunc=lambda x: math.log10(len(x.unique()))) # the log of the count of unique smiles_str
    df.plot.hist(ax=ax,
                 bins = num_bins,
                 title=title,
                 figsize = (15,15))
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    plt.savefig(f"output/{query_results_folder}/{generation_num}/{file_name}_histogram.png")


def plot_scatter(query_results_folder,
                 generation_num,
                 file_name,
                 statistic_col_name,
                 title,
                 x_label,
                 y_label):
    fig, ax = plt.subplots()
    df = pd.read_csv(f"output/{query_results_folder}/{generation_num}/{file_name}.csv")
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
    plt.savefig(f"output/{query_results_folder}/{generation_num}/{file_name}_scatter.png")


def network_statistics(generation_num, query_results_folder):
    """
    Get some statistics on the network.
    0. Number of nodes and edges in the graph, as well as various network-level
        statistics: 1. Eigenvector centrality, 2. Betweenness Centrality,
        3. Random-walk betweenness, 4. Clique enumeration,
        5. k-plex enumeration, 6. k-core enumeration,
        7. k-component enumeration, 8. neighbor redundancy
    1. Node degree distribution: log10 of node degree frequency by degree
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
                                       RETURN algo.asNode(nodeId).smiles_str AS smiles_str, algo.asNode(nodeId).generation_formed AS generation_formed, score AS eigenvector_centrality
                                       ORDER BY eigenvector_centrality DESC """).data()
    save_query_results(generation_num, eigenvector_centrality, "eigenvector_centrality", query_results_folder)
    plot_hist(query_results_folder = query_results_folder,
              generation_num = generation_num,
              file_name = "eigenvector_centrality",
              statistic_col_name = "eigenvector_centrality",
              title = "Histogram of Eigenvector Centrality",
              x_label = "Eigenvector Centrality Score Bin",
              y_label = "Count of Molecules")
    plot_scatter(query_results_folder = query_results_folder,
                 generation_num = generation_num,
                 file_name = "eigenvector_centrality",
                 statistic_col_name = "eigenvector_centrality",
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
                                       RETURN molecule.smiles_str AS smiles_str, molecule.generation_formed AS generation_formed, centrality AS betweenness_centrality
                                       ORDER BY betweenness_centrality DESC;
                                       """).data()
    save_query_results(generation_num, betweenness_centrality, "betweenness_centrality", query_results_folder)
    plot_hist(query_results_folder = query_results_folder,
              generation_num = generation_num,
              file_name = "betweenness_centrality",
              statistic_col_name = "betweenness_centrality",
              title = "Histogram of Betweenness Centrality",
              x_label = "Betweenness Centrality Score Bin",
              y_label = "Count of Molecules")
    plot_scatter(query_results_folder = query_results_folder,
                 generation_num = generation_num,
                 file_name = "betweenness_centrality",
                 statistic_col_name = "betweenness_centrality",
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
                                        RETURN molecule.smiles_str AS smiles_str, molecule.generation_formed AS generation_formed, centrality AS random_walk_betweenness
                                        ORDER BY random_walk_betweenness DESC;""").data()
    save_query_results(generation_num, random_walk_betweenness, "random_walk_betweenness", query_results_folder)
    plot_hist(query_results_folder = query_results_folder,
              generation_num = generation_num,
              file_name = "random_walk_betweenness",
              statistic_col_name = "random_walk_betweenness",
              title = "Histogram of Random Walk Betweenness Centrality",
              x_label = "Random Walk Betweenness Centrality Score Bin",
              y_label = "Count of Molecules")
    plot_scatter(query_results_folder = query_results_folder,
                 generation_num = generation_num,
                 file_name = "random_walk_betweenness",
                 statistic_col_name = "random_walk_betweenness",
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
    RETURN n.smiles_str AS smiles_str, n.exact_mass AS exact_mass,
    n.generation_formed AS generation_formed, size((n)--()) AS count_relationships
    """
    node_deg_query_results = graph.run(node_deg_query).data()
    node_deg_file = "node_distribution_results"
    save_query_results(generation_num = generation_num,
                       query_result = node_deg_query_results,
                       file_name = node_deg_file,
                       this_out_folder = query_results_folder)
    
    # now read in the results, transform, and plot
    # also can represent this as a histogram?
    fig, ax = plt.subplots()
    node_deg_df = pd.read_csv(f"output/{query_results_folder}/{generation_num}/{node_deg_file}.csv")
    # node_deg_df['count_relationships'].value_counts().plot(ax=ax,
    #                                                        kind='bar',
    #                                                        title="Node Degree Distribution by Generation Formed")
    node_deg_pivot = pd.pivot_table(node_deg_df,
                                    values="smiles_str",
                                    index=["count_relationships"],
                                    columns=["generation_formed"],
                                    aggfunc=lambda x: math.log10(len(x.unique()))) # the log of the count of unique smiles_str
    node_deg_pivot.plot(ax=ax,
                        kind="bar",
                        title="Square of Molecule Degree by Generation Formed",
                        figsize = (8,5))
    ax.set_xlabel("Molecule Degree (count of incoming and outgoing edges)")
    ax.set_ylabel("log10(Count of Molecules)")
    plt.savefig(f"output/{query_results_folder}/{generation_num}/{node_deg_file}.png")
    
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
    plt.savefig(f"output/{query_results_folder}/{generation_num}/{node_deg_file}_avg.png")
    
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
    save_query_results(generation_num = generation_num,
                       query_result = incoming_rels_count,
                       file_name = incoming_rels_count_file,
                       this_out_folder = query_results_folder)
    fig, ax = plt.subplots()
    # node_deg_df = pd.read_csv(f"output/{query_results_folder}/{generation_num}/{incoming_rels_count_file}.csv")
    node_deg_df = read_query_results(f"output/{query_results_folder}/{generation_num}/{incoming_rels_count_file}.csv")
    # node_deg_df['count_relationships'].value_counts().plot(ax=ax,
    #                                                        kind='bar',
    #                                                        title="Node Degree Distribution by Generation Formed")
    if not node_deg_df.empty:
        node_deg_pivot = pd.pivot_table(node_deg_df,
                                        values="smiles_str",
                                        index=["count_incoming"],
                                        columns=["generation_formed"],
                                        aggfunc=lambda x: math.log10(len(x.unique()))) # the square of the count of unique smiles_str
        node_deg_pivot.plot(ax=ax,
                            kind="bar",
                            title="Square of Molecule Degree by Generation Formed for Incoming Relationships",
                            figsize = (8,5))
        ax.set_xlabel("Molecule Degree (count of incoming edges)")
        ax.set_ylabel("log10(Count of Molecules)")
        plt.savefig(f"output/{query_results_folder}/{generation_num}/{incoming_rels_count_file}.png")
    
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
    save_query_results(generation_num = generation_num,
                       query_result = outgoing_rels_count,
                       file_name = outgoing_rels_count_file,
                       this_out_folder = query_results_folder)
    fig, ax = plt.subplots()
    # node_deg_df = pd.read_csv(f"output/{query_results_folder}/{generation_num}/{outgoing_rels_count_file}.csv")
    node_deg_df = read_query_results(f"output/{query_results_folder}/{generation_num}/{outgoing_rels_count_file}.csv")
    if not node_deg_df.empty:
        node_deg_pivot = pd.pivot_table(node_deg_df,
                                        values="smiles_str",
                                        index=["count_outgoing"],
                                        columns=["generation_formed"],
                                        aggfunc=lambda x: math.log10(len(x.unique()))) # the square of the count of unique smiles_str
        node_deg_pivot.plot(ax=ax,
                            kind="bar",
                            title="Square of Molecule Degree by Generation Formed for Outgoing Relationships",
                            figsize = (8,5))
        ax.set_xlabel("Molecule Degree (count of outgoing edges)")
        ax.set_ylabel("log10(Count of Molecules)")
        plt.savefig(f"output/{query_results_folder}/{generation_num}/{outgoing_rels_count_file}.png")
    
    # close all plots so they don't accumulate memory
    print("\tNetwork statistics done.")
    plt.close('all')


def graph_from_cypher(data):
    """
    Setting FULL_NETWORK_VISUALIZATION to False because we generated a plot in
    Gephi for the whole network visualizations; not needed in this module. Only
    keeping in case we want to programmatically generate a static network
    visualization.
    
    From: https://stackoverflow.com/questions/59289134/constructing-networkx-graph-from-neo4j-query-result
    
    Constructs a networkx graph from the results of a neo4j cypher query.
    Example of use:
    >>> result = session.run(query)
    >>> G = graph_from_cypher(result.data())

    Nodes have fields 'labels' (frozenset) and 'properties' (dicts). Node IDs correspond to the neo4j graph.
    Edges have fields 'type_' (string) denoting the type of relation, and 'properties' (dict).
    
    """

    G = nx.MultiDiGraph()
    def add_node(node):
        # Adds node id it hasn't already been added
        # print(node)
        # print(type(node))
        # print(node.keys())
        u = node['smiles_str'] # unique identifier for Node
        if G.has_node(u):
            return
        G.add_node(u, labels=node._labels, properties=dict(node))

    def add_edge(relation):
        # Adds edge if it hasn't already been added.
        # Make sure the nodes at both ends are created
        for node in (relation.start_node, relation.end_node):
            add_node(node)
        # Check if edge already exists
        u = relation.start_node['smiles_str'] # unique identifier for Node
        v = relation.end_node['smiles_str'] # unique identifier for Node
        eid = relation['rxn_id'] # unique identifier for Relationship
        if G.has_edge(u, v, key=eid):
            return
        # If not, create it
        G.add_edge(u, v, key=eid, type_=relation.type, properties=dict(relation))

    for d in data:
        for entry in d.values():
            # Parse node
            if isinstance(entry, Node):
                add_node(entry)

            # Parse link
            elif isinstance(entry, Relationship):
                add_edge(entry)
            else:
                raise TypeError("Unrecognized object")
    return G

# def nx2gt(nxG):
#     """
#     Converts a networkx graph to a graph-tool graph.
# 
#     Thanks to: https://bbengfort.github.io/snippets/2016/06/23/graph-tool-from-networkx.html
#     """
#     # Phase 0: Create a directed or undirected graph-tool Graph
#     gtG = gt.Graph(directed=nxG.is_directed())
# 
#     # Add the Graph properties as "internal properties"
#     for key, value in nxG.graph.items():
#         # Convert the value and key into a type for graph-tool
#         tname, value, key = get_prop_type(value, key)
# 
#         prop = gtG.new_graph_property(tname) # Create the PropertyMap
#         gtG.graph_properties[key] = prop     # Set the PropertyMap
#         gtG.graph_properties[key] = value    # Set the actual value
# 
#     # Phase 1: Add the vertex and edge property maps
#     # Go through all nodes and edges and add seen properties
#     # Add the node properties first
#     nprops = set() # cache keys to only add properties once
#     for node, data in nxG.nodes_iter(data=True):
# 
#         # Go through all the properties if not seen and add them.
#         for key, val in data.items():
#             if key in nprops: continue # Skip properties already added
# 
#             # Convert the value and key into a type for graph-tool
#             tname, _, key  = get_prop_type(val, key)
# 
#             prop = gtG.new_vertex_property(tname) # Create the PropertyMap
#             gtG.vertex_properties[key] = prop     # Set the PropertyMap
# 
#             # Add the key to the already seen properties
#             nprops.add(key)
# 
#     # Also add the node id: in NetworkX a node can be any hashable type, but
#     # in graph-tool node are defined as indices. So we capture any strings
#     # in a special PropertyMap called 'id' -- modify as needed!
#     gtG.vertex_properties['id'] = gtG.new_vertex_property('string')
# 
#     # Add the edge properties second
#     eprops = set() # cache keys to only add properties once
#     for src, dst, data in nxG.edges_iter(data=True):
# 
#         # Go through all the edge properties if not seen and add them.
#         for key, val in data.items():
#             if key in eprops: continue # Skip properties already added
# 
#             # Convert the value and key into a type for graph-tool
#             tname, _, key = get_prop_type(val, key)
# 
#             prop = gtG.new_edge_property(tname) # Create the PropertyMap
#             gtG.edge_properties[key] = prop     # Set the PropertyMap
# 
#             # Add the key to the already seen properties
#             eprops.add(key)
# 
#     # Phase 2: Actually add all the nodes and vertices with their properties
#     # Add the nodes
#     vertices = {} # vertex mapping for tracking edges later
#     for node, data in nxG.nodes_iter(data=True):
# 
#         # Create the vertex and annotate for our edges later
#         v = gtG.add_vertex()
#         vertices[node] = v
# 
#         # Set the vertex properties, not forgetting the id property
#         data['id'] = str(node)
#         for key, value in data.items():
#             gtG.vp[key][v] = value # vp is short for vertex_properties
# 
#     # Add the edges
#     for src, dst, data in nxG.edges_iter(data=True):
# 
#         # Look up the vertex structs from our vertices mapping and add edge.
#         e = gtG.add_edge(vertices[src], vertices[dst])
# 
#         # Add the edge properties
#         for key, value in data.items():
#             gtG.ep[key][e] = value # ep is short for edge_properties
# 
#     # Done, finally!
#     return gtG


def network_visualization_by_gen(query_results_folder, generation_num):
    print("Generating an image for the network visualization...")
    # driver = GraphDatabase.driver(url)
    full_network_query = """
    MATCH (n)-[r]->(m)
    RETURN *
    """
    # with driver.session() as session:
    #     result = session.run(full_network_query)
    result = graph.run(full_network_query)
    
    # plot using NetworkX graph object + Matplotlib
    nxG = graph_from_cypher(result.data())
    nx.draw(nxG)
    plt.savefig(f"output/{query_results_folder}/{generation_num}/network_visualization_nxG_at_gen_{generation_num}.png")
    plt.close('all')
    
    # plot using graph_tool module (convert from NetworkX graph to this graph)
    # gtG = nx2gt(nxG)
    # graph_draw(gtG,
    #            vertex_text = g.vertex_index,
    #            output = f"output/{query_results_folder}/{generation_num}/network_visualization_gtG_at_gen_{generation_num}.png")



def compute_likely_abundance_by_molecule(generation_num, query_results_folder):
    """
    Get a dataset with the following columns in order to compute the abundance
    score:
    """
    print("\tComputing the likely abundance score by molecule...")
    # Join all the datasets for rels for all generations. Start with the
    # node_distribution_results query and then join all the other data onto it
    datasets = {'node_distribution_results': ['smiles_str',
                                              'exact_mass',
                                              'generation_formed',
                                              'count_relationships'],
                'incoming_rels_count': ['smiles_str',
                                        'count_incoming'],
                'outgoing_rels_count': ['smiles_str',
                                        'count_outgoing'],
                'betweenness_centrality': ['smiles_str',
                                           'betweenness_centrality'],
                'eigenvector_centrality': ['smiles_str',
                                           'eigenvector_centrality'],
                'random_walk_betweenness': ['smiles_str',
                                            'random_walk_betweenness']
                }
    full_df = pd.DataFrame()
    for dataset in datasets.keys():
        try:
            df = pd.read_csv(f"output/{query_results_folder}/{generation_num}/{dataset}.csv")
            df = df[datasets[dataset]] # filter only by the needed columns
            if dataset == "node_distribution_results":
                full_df = df
            else:
                full_df = pd.merge(full_df, df, on="smiles_str", how='left')
        except:
            # If the df is empty (such as if COLLECT_NETWORK_STATISTICS = False
            # which will make some datasets not exist), then pass
            pass
    
    # Query for all relationships and their parent (consumed) & child (produced)
    # molecules. Groupby molecule and reaction, mulitiply each
    # generation_formed of parent to get likely_abundance_score by relationship.
    # For each molecule, only take most likely reaction (filter out lower scores)
    
    # calculate abundance score
    
    # save calculated data
    full_df.to_csv(f"output/{query_results_folder}/{generation_num}/likely_abundance_score_by_molecule.csv",
                   index=False)
    
    # generate visualization
    
    # save visualization
    pass

# def export_graph_database():
#     """
#     Export nodes / rels to CSV.
#     """
#     pass

def take_network_snapshot(generation_num, query_results_folder, mod_exports_folder_path):
    # do pattern match query on possible autocatalytic cycles
    if PATTERN_MATCHES:
        get_tabulated_possible_autocatalytic_cycles(generation_num = generation_num,
                                                    mod_exports_folder_path = mod_exports_folder_path,
                                                    this_out_folder = query_results_folder,
                                                    ring_size_range = (3, 5),
                                                    feeder_molecule_generation_range = None,
                                                    num_structures_limit = NUM_STRUCTURES_LIMIT) # set to 100 for small batch testing
    
    # do network statistics and get plots
    if COLLECT_NETWORK_STATISTICS:
        network_statistics(generation_num = generation_num,
                           query_results_folder = query_results_folder)
    
    # # get Cytoscape network visualization
    if FULL_NETWORK_VISUALIZATION:
        network_visualization_by_gen(query_results_folder = query_results_folder,
                                     generation_num = generation_num)
    
    # compute likely abundance by molecule
    compute_likely_abundance_by_molecule(generation_num = generation_num,
                                         query_results_folder = query_results_folder)


def get_node_degree_rank_by_gen(df, query_results_folder):
    # Figure #7
    # first rank all molecules by node degree over generation, sorting
    # highest to lowest over generation
    # count_relationships is total count incoming and count outgoing--use this
    # for ranking node degree (lowest rank is highest node degree)
    df['rank_node_deg_by_gen'] = df.groupby(by=['snapshot_generation_num'])['count_relationships'].rank(ascending=False,
                                                                                                        method="max") # so rank values are integer only
    
    
    # plot the results: count_relationships vs rank_node_deg_by_gen coloring
    # by snapshot_generation_num
    # print(df.head())
    # now plot fig 7
    # fig, ax = plt.subplots()
    # df = pd.pivot_table(df,
    #                     values="count_relationships",
    #                     index=["rank_node_deg_by_gen"],
    #                     columns=["snapshot_generation_num"],
    #                     aggfunc=lambda x: math.log10(len(x.unique()))) # the log of the count of unique smiles_str
    # print("pivoted df:")
    # print(df.head())
    # df.plot.area(ax = ax,
    #              stacked = True,
    #              title = "Log of Node Degree by Node Degree Rank by Generation",
    #              figsize = (15,15))
    # ax.set_xlabel("Node Degree Rank")
    # ax.set_ylabel("Log of Node Degree")
    # plt.savefig(f"output/{query_results_folder}/node_deg_by_rank_plot.png")
    
    # df.plot(x = 'rank_node_deg_by_gen',
    #         y = 'count_relationships',
    #         )
    
    # df.plot.area(x='rank_node_deg_by_gen',
    #              y='count_relationships')
    return df


def get_change_in_node_degree_by_molecule_per_gen(df, query_results_folder):
    
    # group by molecule, sort by snapshot generation asc, then take rolling difference
    # of count_relationships
    
    # df.sort_values(by = ['snapshot_generation_num'],
    #                ascending = True,
    #                inplace = True)
    # print(df.head())
    # df['count_rels_diff_by_mol_by_gen'] = df.groupby(by=['smiles_str','snapshot_generation_num'])['count_relationships'].mean().diff()
    # print(df.head())
    
    # calculate the difference in count_relationships by molecule over snapshot generation
    df_mols_with_rels_diff = pd.DataFrame()
    mols = list(df['smiles_str'].unique())
    for mol in mols:
        df_mol = df[ df['smiles_str'] == mol ]
        df_mol.sort_values(by = ['snapshot_generation_num'],
                           ascending = True,
                           inplace = True)
        df_mol['count_rels_diff_by_mol_by_gen'] = df_mol['count_relationships'].diff().fillna(0)
        df_mol = df_mol[['smiles_str', 'snapshot_generation_num', 'count_rels_diff_by_mol_by_gen']]
        df_mols_with_rels_diff = pd.concat([df_mols_with_rels_diff, df_mol])
    
    # Now join this information back onto original df. The calculation
    # count_rels_diff_by_mol_by_gen is computed at the composite level of
    # smiles_str and snapshot_generation_num.
    df = pd.merge(df,
                  df_mols_with_rels_diff,
                  how = "left",
                  on = ['smiles_str', 'snapshot_generation_num'])
    
    # plot the results: count_rels_diff_by_mol_by_gen vs rank_node_deg_by_gen coloring
    # by snapshot_generation_num
    
    
    # finally, return the df
    return df


def get_ring_rels_rule_sequence_motifs(df):
    """
    Example ringRels cell entry:
    [{'rxn_id': '31', 'generation_formed': 1, 'rule': "'Hemiacetal Formation for 7 membered rings'"}, {'rxn_id': '254', 'rule': "'Elimination + enol to keto'", 'generation_formed': 1}, {'rxn_id': '126', 'rule': "'Cannizarro 2", 'generation_formed': 1}, {'rxn_id': '155', 'rule': "'Knoevenagel H (inv)'", 'generation_formed': 1}, {'rxn_id': '40', 'rule': "'Aldol Condensation'", 'generation_formed': 1}]
    """
    # convert string representation of list of dicts to objects
    # all_rings = list(df['ringRels'].unique())
    # all_rings_cleaned = []
    # for ring_str in all_rings:
    #     # cleaned_str = ring_str.replace("'","\"\"")
    #     ring_rels = json.loads(ring_str)
    #     all_rings_cleaned.append(ring_rels)
    # 
    # # record all rule sequences and count of how many times the sequence appears
    # rule_sequences = {} # hold counts of rule sequences
    # for ring in all_rings_cleaned:
    #     for rel in ring:
    #         print(rel['rule'])
    #     wait = input("Press enter...")
    pass
    


def compile_all_generations_data(query_results_folder, generation_limit):
    print("\t\tCompiling abundance score data and autocatalysis pattern query results from all network snapshots into one file...")
    
    # pull all molecule scores calculated at each generation thus far
    out_dir = f"output/{query_results_folder}"
    df_all_gens = pd.DataFrame()
    df_autocat_all_gens = pd.DataFrame()
    for generation_num in range(generation_limit + 1):
        # append molecule data
        gen_file_path = out_dir + f"/{generation_num}/likely_abundance_score_by_molecule.csv"
        try:
            df_gen = pd.read_csv(gen_file_path)
            df_gen['snapshot_generation_num'] = generation_num
            df_all_gens = pd.concat([df_all_gens, df_gen])
        except:
            pass
        
        # append autocat query results data
        autocat_gen_file_path = out_dir + f"/{generation_num}/autocat_query_results.csv"
        try:
            df_autocat_gen = pd.read_csv(autocat_gen_file_path)
            df_autocat_gen['snapshot_generation_num'] = generation_num
            df_autocat_all_gens = pd.concat([df_autocat_all_gens, df_autocat_gen])
        except:
            # if failed, pd.read_csv was an empty file (no autocat pattern matches found)
            pass
    
    # do some more calculations/analysis/plotting on all generations' dataset
    df_all_gens = get_node_degree_rank_by_gen(df_all_gens, query_results_folder)
    df_all_gens = get_change_in_node_degree_by_molecule_per_gen(df_all_gens, query_results_folder)
    
    # save all generations' molecule data as CSV
    df_all_gens.to_csv(out_dir + "/all_generations_abundance_scores.csv",
                       index=False)
    
    # save all generations' autocatalysis pattern matches as CSV
    df_autocat_all_gens.to_csv(out_dir + "/all_generations_autocat_pattern_matches.csv",
                               index = False)
    
    # get rule sequence motifs (frequently found sequences in the ring relationships)
    get_ring_rels_rule_sequence_motifs(df_autocat_all_gens)
    
    
def smiles_passes_filter(smiles_str):
    passes_filter = True
    for filt_mol in MOLECULE_FILTER:
        if filt_mol == smiles_str:
            passes_filter = False
    return passes_filter
    

def split_rels_into_from_and_to_mols(rels):
    """
    To parse Aayush's and Romulo's rels output format.
    """
    from_mols = []
    to_mols = []
    for rel in rels:
        if rel != '':
            rel_data = rel.split('\t')
            gen_consumpt_scalar = int(rel_data[2])
            if gen_consumpt_scalar == -1:
                # generation/consumption scalar is negative; this species must be
                # being consumed in the reaction
                from_mols.append(rel)
            elif gen_consumpt_scalar == 1:
                # generation/consumption scalar is positive; this species must be
                # being generated in the reaction
                to_mols.append(rel)
    return from_mols, to_mols

def filter_to_mols_by_rxn_id(rxn_id, to_mols):
    filtered_to_mols = []
    for rel in to_mols:
        rel_data = rel.split("\t")
        this_rxn_id = rel_data[0]
        if this_rxn_id == rxn_id:
            filtered_to_mols.append(rel)
    return filtered_to_mols



def import_data_from_MOD_exports(mod_exports_folder_path, network_name, generation_limit):
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
    nodes_all_generation_file_names = os.listdir(nodes_folder)
    rels_all_generation_file_names = os.listdir(rels_folder)
    
    # first get all generation numbers from the Neo4j_Imports .txt files
    # and sort them in order so the import doesn't go out of order
    all_gens = []
    for generation_file in nodes_all_generation_file_names:
        generation_num = int(generation_file.split("_")[1].split(".")[0])
        all_gens.append(generation_num)
    all_gens.sort()
    max_all_gens = max(all_gens)
    # all_gens.append(max_all_gens + 1) # for rels_import_gen shift; explained below
    
    # if no given generation_limit, set the generation limit to the
    # max (so no data type issue between None/Integer)
    if generation_limit == None:
        generation_limit = max_all_gens
    # don't let the generation limit go above the max exported from MOD
    if generation_limit > max_all_gens:
        generation_limit = max_all_gens
    
    # set up the output folder for snapshots (statistics/visualizations) of the network
    query_results_folder = network_name + "_" + get_timestamp()
    # query_results_folder = "2020-08-01_18-39-23-986232" # manually override folder name for debugging
    os.mkdir("output/" + query_results_folder)
    # Optional: save the state of the Neo4j_Imports folder at the time this was run
    copytree(mod_exports_folder_path, 'output/' + query_results_folder + "/Neo4j_Imports")
    
    # now import the data.
    for generation_num in all_gens:
        # Shift the rels data import number by 1 because the nodes_1.txt file
        # actually represents the 0th generation. So when nodes_1.txt file imports,
        # no rels file will import, then when nodes_2.txt file imports, rels_1.txt
        # will import, and so on.
        # rels_import_gen = generation_num - 1
        if generation_num <= generation_limit:
            print(f"\tGeneration number {generation_num}...")
            # load in this generation's data
            rels_generation_file = "rels_" + str(generation_num) + ".txt"
            rels = open(rels_folder + "/" + rels_generation_file,'r').read().split('\n')
            # if SHUFFLE_GENERATION_DATA is implemented, shuffle the order
            # of the molecules to be imported
            if SHUFFLE_GENERATION_DATA:
                random.shuffle(rels) # random.shuffle(nodes)
            # from_mols, to_mols = split_rels_into_from_and_to_mols(rels)
            
            # create an output folder for this generation within query_results_folder
            os.mkdir("output/" + query_results_folder + f"/{generation_num}")
            
            # import molecule nodes
            print("\t\tImporting Molecule and Reaction nodes...")
            
            # iterate through molecules (make sure all molecule nodes are imported
            # before we try to merge edges onto them)
            for rel in rels:
                if rel != '':
                    rel_data = rel.split('\t')
                    # import molecule node
                    smiles_str = rel_data[1]
                    if smiles_str != "":
                        create_molecule_if_not_exists(smiles_str = smiles_str,
                                                      generation_formed = generation_num)
                    # import reaction node
                    rxn_id = rel_data[0]
                    rxn_rule = rel_data[3]
                    if rxn_id != "":
                        create_reaction_if_not_exists(id = rxn_id,
                                                      rule = rxn_rule,
                                                      generation_formed = generation_num)
            # wait = input("Press enter...")
            
            # create relationship edges
            print("\t\tImporting REACTANT/PRODUCT edges...")
            for rel in rels:
                if rel != '':
                    rel_data = rel.split('\t')
                    rxn_id = rel_data[0]
                    smiles_str = rel_data[1]
                    gen_consumpt_scalar = int(rel_data[2])
                    rxn_rule = rel_data[3]
                    if gen_consumpt_scalar == 1:
                        # this molecule is being generated, therefore it is a PRODUCT
                        create_product_rel_if_not_exists(smiles_str = smiles_str,
                                                         rxn_id = rxn_id,
                                                         generation_formed = generation_num)
                    elif gen_consumpt_scalar == -1:
                        # this molecule is being consumed, therefore it is a REACTANT
                        create_reactant_rel_if_not_exists(smiles_str = smiles_str,
                                                          rxn_id = rxn_id,
                                                          generation_formed = generation_num)
            # wait = input("Press enter...")
            
            # Now that the generation's data has been loaded into the network,
            # take a snapshot of it. Only take snapshot at each generation,
            # if NETWORK_SNAPSHOTS is True. Otherwise, only take a snapshot
            # if this is the last generation.
            if ((not NETWORK_SNAPSHOTS) and (generation_num == generation_limit)):
                print("\t\tTaking statistics/visualizations snapshot of network...")
                take_network_snapshot(generation_num,
                                      query_results_folder,
                                      mod_exports_folder_path)
            elif NETWORK_SNAPSHOTS:
                print("\t\tTaking statistics/visualizations snapshot of network...")
                take_network_snapshot(generation_num,
                                      query_results_folder,
                                      mod_exports_folder_path)
    
    # now that all snapshots of the network have been taken, compile all of the
    # data into one source (only if NETWORK_SNAPSHOTS enabled)
    # if NETWORK_SNAPSHOTS:
    #     compile_all_generations_data(query_results_folder,
    #                                  generation_limit)


if __name__ == "__main__":
    for _ in range(REPEAT_RUNS):
        for mod_export_folder_path in EXPORT_PATHS:
            print(f"Importing the network from the following path: {mod_export_folder_path}")
            network_name = mod_export_folder_path.split('/')[-2]
            import_data_from_MOD_exports(mod_exports_folder_path = mod_export_folder_path,
                                         network_name = network_name,
                                         generation_limit = GENERATION_LIMIT) # Set to None or Integer. The generation limit at which to import
    
    # test functions
    # query_results_folder = "2020-10-15_17-35-22-778323"
    # generation_num = 2
    # network_visualization_by_gen(query_results_folder, generation_num)












