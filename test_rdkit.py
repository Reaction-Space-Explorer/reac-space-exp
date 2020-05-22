# -*- coding: utf-8 -*-
from rdkit import Chem
from rdkit.Chem import AllChem
import networkx as nx
import matplotlib.pyplot as plt
import random
from models import *
import pandas as pd
import itertools
import datetime
# import requests
import re
import codecs

# enable this bit if you want lots of print output
debug = False

def timestamp():
    """
    Create timestamp safe for use in file names.
    """
    illegal_chars = "* . \" / \\ [ ] : ; | ,".split(" ")
    timestamp_str = datetime.datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
    for char_str in illegal_chars:
        timestamp_str = timestamp_str.replace(char_str, "-")
    return timestamp_str


def clear_db_table():
    """
    Delete all records from table and start over each run.
    If keeping data for multiple runs, just add column run_id which
    increments at the start of every run.
    """
    execute_query(sql="DELETE FROM molecule;", params=())
    

def react(reactant_a, reactant_b):
    """
    Receive list of reactants (list size 2), return list of products (variable size).
    
    RDKit Chemical Reactions documentation:
    https://www.rdkit.org/docs/GettingStartedInPython.html#chemical-reactions
    
    SMARTS examples:
    https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
    
    Need to figure out list SMARTS formats for ReactionFromSmarts, assuming that
    multiple are needed. If only one is needed, can remove for loop. Could also
    take this list as an input / set as global var if needs to change.
    """
    smarts_format_rxns = ['[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]',
                          '[C:1]=[O,N:2]>>[C:1][*:2]',
                          '[C:1]=[O,N:2]>>*[C:1][*:2]',
                          '[C:1]~[O,N:2]>>*[C:1]~[*:2]',
                          "([C:1]=[C;H2].[C:2]=[C;H2])>>[*:1]=[*:2]",
                          '[CH1:1][OH:2].[OH][C:3]=[O:4]>>[C:1][O:2][C:3]=[O:4]',
                          '[C@H1:1][OH:2].[OH][C:3]=[O:4]>>[C@:1][O:2][C:3]=[O:4]',
                          '(C(=O)O).(OCC)>>C(=O)OCC.O',
                          '[C:1]=[C:2] + [C:3]=[C:4][C:5]=[C:6] -> [C:1]1[C:2][C:3][C:4]=[C:5][C:6]1']
    all_products = []
    for smart_str in smarts_format_rxns:
        try:
            # try to react these reactants
            # print(reactant_a)
            # print(reactant_b)
            rxn = AllChem.ReactionFromSmarts(smart_str)
            # RunReactants accepts a tuple, so this might be any length (more than 2 reactants possible?)
            # default maxProducts for RunReactants is 1000
            product_set = rxn.RunReactants(reactants=(reactant_a, reactant_b))
            # print("PRODUCT_SET: ", product_set)
            # for product in product_set:
            #     print('PRODUCT: ', product)
            product_set = [product_set[0][0] for prod in product_set] # for some reason RunReactants returns a tuple of tuple size 2, where only the first index contains a product...
            # print("PRODUCT_SET: ", product_set)
            all_products.extend(product_set) # add list of products from each reaction to big all_products list
        except:
            # the products in this reaction will not react
            # or there is an error with the smarts_str
            pass
    all_products = list(set(all_products)) # remove possible duplicate products with set(), then convert back to list
    # print("ALL PRODUCTS:", all_products)
    return all_products


def reaction_already_computed(reactant_a, reactant_b):
    """
    No need to duplicate reaction calculations.
    
    Check molecule table in network.db to determine whether or not reaction has
    already been computed in a prior generation.
    """
    s = Session()
    # look up ids of reactants
    if (reactant_a != None) and (reactant_b != None):
        r_a = s.query(Molecule).filter(Molecule.smiles_str == Chem.MolToSmiles(reactant_a)).first()
        r_b = s.query(Molecule).filter(Molecule.smiles_str == Chem.MolToSmiles(reactant_b)).first()
        
        # two possible orders (reactant_1 in db can either really be reactant_a or _b;
        # order doesn't matter )
        rxn_order_1 = s.query(Molecule).filter(Molecule.reactant_1 == r_a.rowid).filter(Molecule.reactant_2 == r_b.rowid).first()
        rxn_order_2 = s.query(Molecule).filter(Molecule.reactant_1 == r_b.rowid).filter(Molecule.reactant_2 == r_a.rowid).first()
        s.close()
        
        # if both reaction orders are None, reaction has not been computed before
        if (rxn_order_1 == None) and (rxn_order_2 == None):
            return False # reaction not already computed
        else:
            True # reaction already computed, proceed to calculate
    else:
        return None

def create_molecule_record(molecule_obj, generation_count, reactant_1, reactant_2):
    """
    Create new record for a molecule in the db with generation information.
    (Any information available in RDKit molecule object can be stored at this point
    as well; for now all I'm storing is the molecule's unique identifier to pull 
    information about it later.)
    """
    s = Session()
    # print(reactant_1, reactant_2, molecule_obj)
    try:
        r1_m = Chem.MolToSmiles(reactant_1)
    except:
        r1_m = None
        
    try:
        r2_m = Chem.MolToSmiles(reactant_2)
    except:
        r2_m = None
    
    if (reactant_1 != None) and (r1_m != None):
        reactant_1 = s.query(Molecule).filter(Molecule.smiles_str == r1_m).first().rowid
    if (reactant_2 != None) and (r2_m != None):
        reactant_2 = s.query(Molecule).filter(Molecule.smiles_str == r2_m).first().rowid
    new_m = Molecule(reactant_1 = reactant_1,
                     reactant_2 = reactant_2,
                     smiles_str = Chem.MolToSmiles(molecule_obj),
                     generation_formed = generation_count)
    s.add(new_m)
    s.commit()
    s.close()

def molecule_in_db(molecule_obj):
    """
    Check whether molecule has record in network.db table "molecule".
    """
    s = Session()
    # print("molecule_obj: ", molecule_obj)
    # try:
    if molecule_obj != None:
        m = s.query(Molecule).filter(Molecule.smiles_str == Chem.MolToSmiles(molecule_obj)).first()
        if m == None:
            return False
        else:
            return True
    else:
        return None
    # except Exception as e:
    #     print("molecule_in_db: Error with the following molecule_obj:", molecule_obj)
    #     s.close()
    #     return None
    


def calculate_generations(all_reactants, num_generations, generation_count):
    print(f"\t{generation_count}...")
    len_reactants_start = len(all_reactants)
    # if generation 0, then no molecules in db; create records for initial set without reactants
    if generation_count == 0:
        if debug: print("\t\tRecording input molecules...")
        for m in all_reactants:
            exists = molecule_in_db(m)
            if (exists != None) and (not exists): # None means molecule_in_db threw an error
                create_molecule_record(molecule_obj = m,
                                       generation_count = generation_count,
                                       reactant_1 = None,
                                       reactant_2 = None)
    
    # compute list of all possible reaction combinations from global reactants list
    if debug: print("\t\tFinding all possible reactant pair combinations...")
    reactant_combinations = list(itertools.combinations(all_reactants, 2)) # assume always 2 reactants per reaction
    if debug: print("ALL REACTANTS:", all_reactants)
    if debug: wait = input('Press enter to continue or Ctrl+C to exit...')
    if debug: print("ALL REACTANT COMBINATIONS:",reactant_combinations)
    if debug: wait = input('Press enter to continue or Ctrl+C to exit...')
    
    # loop through reactant combinations to call react() and get product data
    if debug: print("\t\tComputing any new reactions and logging information on any new products formed...")
    num_prods_this_gen = 0
    for rxn in reactant_combinations:
        reactant_a = rxn[0]
        reactant_b = rxn[1]
        if not reaction_already_computed(reactant_a, reactant_b):
            products = react(reactant_a = reactant_a, reactant_b = reactant_b)
            if debug: print(products)
            for product in products:
                if not molecule_in_db(product):
                    if debug: print(product) # print out only new products
                    # record generation at which product was formed, as well as reactants
                    create_molecule_record(molecule_obj = product,
                                           generation_count = generation_count,
                                           reactant_1 = reactant_a,
                                           reactant_2 = reactant_b)
                if product not in all_reactants:
                    all_reactants.append(product) # insert product as reactant for next generation
                    num_prods_this_gen += 1 # for break condition (where if no new products are formed this generation, trigger end of script)
        # print(f"\t\tnum_prods_this_gen: {num_prods_this_gen}")
            # rxn_data = [reactant_a, reactant_b, product] # insert other data of interest from a reaction here in a data dict {}
            # rxn_outcomes.append(rxn_data) # append data on newly processed reactions
            # all_rxns.append(rxn_data[:2]) # append meta data only onto all_rxns list (right now no reaction data associated, so only reaction + product species names available)
    if debug: wait = input('Press enter to continue or Ctrl+C to exit...')
    len_reactants_after_prods_computed = len(all_reactants)
    print("All reactants len() before and after computing products:", len_reactants_start, len_reactants_after_prods_computed)
    # if break condition reached, stop, otherwise run next generation
    if debug: print("\t\tChecking if break condition has been reached...")
    generation_count += 1
    if generation_count >= num_generations:
        return print("All done! Reached maximum number of generations to compute.")
    elif num_prods_this_gen == 0:
        return print("All done! Ran out of reactions to do.")
    else:
        calculate_generations(all_reactants = all_reactants,
                              num_generations = num_generations, # stop condition never changes from main(), just keep passing to next generation
                              generation_count = generation_count)

def export_network_to_csv_with_self_joins():
    """
    Export "molecule" table in DB to network.csv with timestamp.
    
    Do self-referential join to join parent (i.e. reactant) information
    onto each product molecule record. (OUTER JOIN because input
    molecule records will not have a parent record for reactant_1 and _2.)
    """
    df = pd.read_sql("""SELECT molecule.reactant_1 AS reactant_1_id,
                     r1.smiles_str AS reactant_1_smiles_str,
                     molecule.reactant_2 AS reactant_2_id,
                     r2.smiles_str AS reactant_2_smiles_str,
                     molecule.smiles_str,
                     molecule.generation_formed
                     FROM molecule
                     LEFT OUTER JOIN molecule AS r1 ON r1.rowid = molecule.reactant_1
                     LEFT OUTER JOIN molecule AS r2 ON r2.rowid = molecule.reactant_2
                     ;""", con=engine)
    df.to_csv("output/network_" + timestamp() + ".csv", index=False)

def export_network_to_csv():
    """
    Export "molecule" table in DB to network.csv with timestamp.
    
    Do self-referential join to join parent (i.e. reactant) information
    onto each product molecule record. (OUTER JOIN because input
    molecule records will not have a parent record for reactant_1 and _2.)
    """
    df = pd.read_sql("""SELECT * FROM molecule;""",
                     con=engine)
    df.to_csv("output/network_" + timestamp() + ".csv", index=False)

def load_rxn_network_from_db():
    """
    Take Molecule table from db and load into graph management system of choice
    (i.e. NetworkX or Neo4j).
    """
    # visualize reaction network
    # nx.draw(rxn_network, with_labels=True, arrows=True)
    # plt.show()
    pass



def scrape_molecules():
    """
    Scrape RDKit documentation for molecules because I'm too lazy to copy & paste.
    """
    # Using requests library for only one GET request, I received the response "Not Acceptable!"... Sorry RDKit...
    # url = "http://www.rdkit.org/docs/RDKit_Book.html#reaction-smarts"
    # txt = requests.get(url).text
    # print(txt[0:30])
    txt = codecs.open("RDKit_book/book_txt.txt",'r','utf-8').read()
    match_str = r"MolFromSmiles\('.*'\)"
    all_mols = re.findall(match_str, txt) # returns a list of strings, such as the string "MolFromSmiles('C=O')"
    # print(len(all_mols))
    all_mols_file = open('input_molecules.txt','w')
    for mol in all_mols:
        mol = mol.split(r"'")[1]
        # print(mol)
        all_mols_file.write(mol + "\n")
    all_mols_file.close()
        


def main():
    """
    Configure and run it!
    """
    print("Loading input molecules...")
    # define input molecules (SMILES format)
    # test mX molecules taken from here: http://www.rdkit.org/docs/RDKit_Book.html#reaction-smarts
    # molecule_def = {"glycine": "C(C(=O)O)N",
    #                 "HCN": "C#N",
    #                 "H2CO": "C=O",
    #                 "glucose": "C(C1C(C(C(C(O1)O)O)O)O)O",
    #                 'm1': 'C1=CC2=C(C=C1)C1=CC=CC=C21',
    #                 "m2": 'O=C1C=CC(=O)C2=C1OC=CO2',
    #                 'm3': 'C1=C[N]C=C1',
    #                 'm4': 'C1=CC=CC=C[C+]1',
    #                 'm5': 'C1=[C]NC=C1',
    #                 'm6': 'OC(=O)c1[te]ccc1',
    #                 'm7': 'CC=O',
    #                 'm8': 'CC=N',
    #                 'm9': 'CC=O',
    #                 'm10': 'C=CCOCC=C',
    #                 'm11': 'CC(CCN)O'
    #                 }
    # # get list of rdkit molecule objects from molecule definitions
    # input_molecules = [Chem.MolFromSmiles(molecule_def[key]) for key in molecule_def.keys()]
    
    # extend hardcoded input list with a bunch of molecules from RDKit documentation
    # scrape_molecules()
    scraped_molecules = codecs.open("input_molecules.txt",'r','utf-8').read().split("\n")
    scraped_molecules = [Chem.MolFromSmiles(mol) for mol in scraped_molecules]
    # input_molecules.extend(scraped_molecules)
    input_molecules = scraped_molecules
    # print(molecules[0]) # should return __repr__ of molecule object; if loading molecule failed, returns None
    
    # clear all records from db table (generated from last time main() was run)
    print("Preparing database...")
    clear_db_table()
    
    # compute all generations, up to num_generations (or until products list is empty,
    # i.e. no new reactions.)
    # initialize generation_count as 0; variable will pass through recursively
    # to increment by 1 each generation
    print("Calculating generation number:")
    try:
        calculate_generations(all_reactants = input_molecules,
                              num_generations = 5,
                              generation_count = 0)
        export_network_to_csv()
    except KeyboardInterrupt:
        export_network_to_csv()
        print("Keyboard Interrupt received. Output exported.")
    # export_network_to_csv_with_self_join() # very slow; just load into network
    
    


# call the function!
main()
# scrape_molecules()






