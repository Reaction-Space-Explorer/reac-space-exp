import pandas as pd
import numpy as np
import os
from py2neo import Graph, Node, Relationship, NodeMatcher, RelationshipMatcher

url = "bolt://neo4j:0000@localhost:7687"
graph = Graph(url)
matcher = NodeMatcher(graph)
rel_matcher = RelationshipMatcher(graph)


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
        


def import_data_from_MOD_exports(folder_path):
    """
    Import simple fake dataset to test queries out on.
    """
    # first, clear db
    print("Clearing database...")
    graph.run("MATCH (n) DETACH DELETE n")
    
    # read in data
    nodes_folder = folder_path + "/nodes"
    rels_folder = folder_path + "/rels"

    # create nodes by iterating through each generation's file
    print("Importing nodes...")
    nodes_all_generation_file_names = os.listdir(nodes_folder)
    for generation_file in nodes_all_generation_file_names:
        generation_num = int(generation_file.split("_")[1].split(".")[0])
        print(f"\tGeneration number {generation_num}...")
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
    for generation_file in rels_all_generation_file_names:
        generation_num = int(generation_file.split("_")[1].split(".")[0])
        print(f"\tGeneration number {generation_num}...")
        rels = open(rels_folder + "/" + generation_file,'r').read().split('\n')
        for rel in rels:
            if rel != "":
                rel_data = rel.split(',')
                create_relationship_if_not_exists(rxn_id = rel_data[0],
                                                  from_smiles = rel_data[1],
                                                  to_smiles = rel_data[2],
                                                  rule = rel_data[3],
                                                  generation_formed = generation_num)
        # merge_query = rxn_query_str(reactant=row['from_node'],
        #                             product=row['to_node'],
        #                             rxn_id=row['rxn_id'])


# load_graph()

# import_molecules()

# import_mock_data()

# mod_exports_folder_path = "../radicals/Neo4j_Imports"
mod_exports_folder_path = "mock_data/Neo4j_Imports_Test"
import_data_from_MOD_exports(mod_exports_folder_path)












