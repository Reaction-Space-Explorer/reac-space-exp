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


# configure network database Neo4j
url = "bolt://neo4j:0000@localhost:7687"
graph = Graph(url)
matcher = NodeMatcher(graph)
rel_matcher = RelationshipMatcher(graph)



def rxn_query_str(reactant, product, rxn_id):
    """
    Generate cypher MERGE query for reactant and product node.
    """
    return "MERGE (r: Molecule {smiles_str:\""+ reactant +"\"}) MERGE (p: Molecule {smiles_str:\""+ product +"\"}) MERGE (r)-[:FORMS {rxn_id: "+ str(rxn_id) +"}]->(p)"


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

def create_reaction_rel_if_not_exists(from_smiles, to_smiles):
    from_molecule = matcher.match("Molecule", smiles_str = from_smiles).first()
    to_molecule = matcher.match("Molecule", smiles_str = to_smiles).first()
    if len(list(graph.match(nodes=(from_molecule, to_molecule), r_type="FORMS"))) <= 0:
        # relationship does not exist
        tx = graph.begin()
        new_r = Relationship(from_molecule, "FORMS", to_molecule)
        tx.create(new_r)
        tx.commit()


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
        create_reaction_rel_if_not_exists(from_smiles = row['from_node'],
                                          to_smiles = row['to_node'])
        # merge_query = rxn_query_str(reactant=row['from_node'],
        #                             product=row['to_node'],
        #                             rxn_id=row['rxn_id'])


# def import_molecules():
#     """
#     Takes all output .txt files from data folder and parses the text files for
#     any new molecules each generation.
#     """
#     txt_file_names = os.listdir(os.path.join(os.getcwd(), "data"))
#     for file_name in txt_file_names:
#         generation = int(file_name.split(".")[0][-1])
#         molecules = open(f"data/molecules/{file_name}").read().split("\n")
#         for molecule in molecules:
#             create_molecule_if_not_exists(smiles_str = molecule,
#                                           generation_formed = generation)


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



# import_mock_data()