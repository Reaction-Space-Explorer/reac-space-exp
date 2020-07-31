import pandas as pd
import numpy as np
import os
from py2neo import Graph, Node, Relationship, NodeMatcher, RelationshipMatcher
import json
import datetime
import matplotlib.pyplot as plt
import shutil

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
                                properties = {"rule": rule, "rxn_id": rxn_id, "generation_formed": generation_formed})
    if len(list(match_pattern)) <= 0:
        # relationship does not exist
        tx = graph.begin()
        # see documentation for weird Relationship function; order of args go:
        # from node, relationship, to node, and then kwargs for relationship properties
        # https://py2neo.org/v4/data.html#py2neo.data.Relationship
        new_r = Relationship(from_molecule, "FORMS", to_molecule, rule=rule, rxn_id=rxn_id, generation_formed = generation_formed)
        tx.create(new_r)
        tx.commit()


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
                     exact_mass = exact_mass,
                     generation_formed = generation_formed)
        tx.create(new_m)
        tx.commit()
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
        


def import_data_from_MOD_exports(mod_exports_folder_path):
    """
    Import simple fake dataset to test queries out on.
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
    for generation_num in all_gens:
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



def get_tabulated_possible_autocatalytic_cycles(mod_exports_folder_path,
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
    # print("\t\t\t" + query_txt)
    
    # execute query in Neo4j
    print("\t\tExecuting query and collecting results (this may take awhile)...")
    query_result = graph.run(query_txt).data()
    # print("\t\tQuery results:")
    # print(query_result[0])
    print("\t\tSaving query results and meta info...")
    this_out_folder = get_timestamp()
    os.mkdir("output/" + this_out_folder)
    # save data as JSON and CSV (JSON for code readability, CSV for human readability)
    with open('output/' + this_out_folder + '/query_results.json', 'w') as file_data_out:
        json.dump(query_result, file_data_out)
    data_df = pd.read_json('output/' + this_out_folder + '/query_results.json')
    data_df.to_csv('output/' + this_out_folder + '/query_results.csv', index=False)
    
    # save Neo4j_Imports folder
    shutil.copyfile(mod_exports_folder_path, 'output/Neo4j_Imports')
    
    # save meta info as well in out folder
    with open("output/" + this_out_folder + "/query.txt", 'w') as file_query_out:
        file_query_out.write(query_txt)
    query_params = pd.DataFrame( {"parameter": ["min_ring_size","max_ring_size","min_feeder_gen","max_feeder_gen","num_structures_limit"],
                                  "value": [min_ring_size, max_ring_size, min_feeder_gen, max_feeder_gen, num_structures_limit] } )
    query_params.to_csv("output/" + this_out_folder + "/query_parameters.csv", index=False)
    return this_out_folder



def analyze_possible_autocatalytic_cycles(mod_exports_folder_path):
    """
    Now that we have the tabulated results of the graph queries, let's do some
    analysis on what's going on.
    
    1. Calculate the the count of cycles found per generation
    2. Total/avg mass per cycle per generation
    3. 
    """
    
    print("Preparing graph query:")
    query_results_folder = get_tabulated_possible_autocatalytic_cycles(mod_exports_folder_path = mod_exports_folder_path,
                                                                       ring_size_range = (3, 6),
                                                                       feeder_molecule_generation_range = None,
                                                                       num_structures_limit = 10)
    # query_results_folder = "2020-07-27_15-48-35-964434" # manually override for debugging
    
    print("Generating some plots on cycle size distribution / stats by generation...")
    query_data = pd.read_json("output/" + query_results_folder + "/query_results.json")
    # print(query_data.describe())
    # print(query_data.head())
    # cycle distribution (y axis is frequency, x axis is ring size)
    fig, ax = plt.subplots()
    query_data['countMolsInRing'].value_counts().plot(ax = ax,
                                                      kind='bar',
                                                      title = "Ring Size Frequency Distribution")
    
    plt.show()
    
    print("Saving plots...")
    
    print("Network analysis done.")






if __name__ == "__main__":
    # mod_exports_folder_path = "../main/Neo4j_Imports"
    mod_exports_folder_path = "../radicals/all7/Neo4j_Imports"
    # import_data_from_MOD_exports(mod_exports_folder_path = mod_exports_folder_path)
    analyze_possible_autocatalytic_cycles(mod_exports_folder_path = mod_exports_folder_path)
    







