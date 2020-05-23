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


def log_reaction(reactant_1, reactant_2, reaction_id, produced_products):
    # get molecule objects for their rowid
    s = Session()
    try:
        m1_id = s.query(Molecule.rowid).filter(Molecule.smiles_str == Chem.MolToSmiles(reactant_1)).first().rowid
    except:
        m1_id = None
    try:
        m2_id = s.query(Molecule.rowid).filter(Molecule.smiles_str == Chem.MolToSmiles(reactant_2)).first().rowid
    except:
        m2_id = None
    
    # create record in db ReactionLog table
    new_rxn_log = ReactionLog(reaction_id = reaction_id,
                              reactant_1 = m1_id,
                              reactant_2 = m2_id,
                              produced_products = produced_products)
    s.add(new_rxn_log)
    s.commit()
    s.close()


def react(reactant_a, reactant_b, rxn_obj, rxn_id):
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
    # log that this reaction was attempted no matter whether successful or not;
    # don't want to try it again
    
    try:
        # try to react these reactants
        # print(reactant_a)
        # print(reactant_b)
        # RunReactants accepts a tuple, so this might be any length (more than 2 reactants possible?)
        # default maxProducts for RunReactants is 1000
        product_set = rxn_obj.RunReactants(reactants=(reactant_a, reactant_b))
        
        
        # for some reason RunReactants returns a tuple of tuple size 2, where only the first index contains a product...
        # set() gets unique list of products, list() converts set object back to list object
        products = list(set([product_set[0][0] for prod in product_set]))
        log_reaction(reactant_1 = reactant_a,
                     reactant_2 = reactant_b,
                     reaction_id = rxn_id,
                     produced_products = 1)
        return products
    except:
        # the products in this reaction will not react
        # or there is an error with the smarts_str of the rxn_obj
        log_reaction(reactant_1 = reactant_a,
                     reactant_2 = reactant_b,
                     reaction_id = rxn_id,
                     produced_products = 0)


def reaction_already_computed(reactant_a, reactant_b, rxn_id):
    """
    No need to duplicate reaction calculations.
    
    Check molecule table in network.db to determine whether or not reaction has
    already been computed in a prior generation.
    """
    s = Session()
    # look up ids of reactants
    if (reactant_a != None) and (reactant_b != None):
        # get molecule database objects
        r_a = s.query(Molecule).filter(Molecule.smiles_str == Chem.MolToSmiles(reactant_a)).first()
        r_b = s.query(Molecule).filter(Molecule.smiles_str == Chem.MolToSmiles(reactant_b)).first()
        
        # two possible orders (reactant_1 in db can either really be reactant_a or _b;
        # order doesn't matter )
        rxn_order_1 = s.query(ReactionLog).filter(ReactionLog.reactant_1 == r_a.rowid).filter(ReactionLog.reactant_2 == r_b.rowid).filter(ReactionLog.reaction_id == rxn_id).first()
        rxn_order_2 = s.query(ReactionLog).filter(ReactionLog.reactant_1 == r_b.rowid).filter(ReactionLog.reactant_2 == r_a.rowid).filter(ReactionLog.reaction_id == rxn_id).first()
        s.close()
        
        # if both reaction orders are None, reaction has not been computed before
        if (rxn_order_1 == None) and (rxn_order_2 == None):
            return False # reaction not already computed, proceed to calculate
        else:
            True # reaction already computed
    else:
        return None # reactant_a or reactant_b cannot be converted to a molecule object; error

def get_or_create_reaction_record(smarts_str):
    s = Session()
    check_exists = s.query(Reaction).filter(Reaction.smarts_str == smarts_str).first()
    if check_exists == None:
        new_r = Reaction(smarts_str = smarts_str)
        s.add(new_r)
        s.commit()
        rowid = new_r.rowid
        s.close()
        return rowid
    else:
        rowid = check_exists.rowid
        return rowid
    



def create_molecule_record(molecule_obj, generation_count, reactant_1, reactant_2):
    """
    Create new record for a molecule in the db with generation information.
    (Any information available in RDKit molecule object can be stored at this point
    as well; for now all I'm storing is the molecule's unique identifier to pull 
    information about it later.)
    """
    if molecule_obj != None:
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
                         generation_formed = generation_count,
                         exact_mass = Chem.rdMolDescriptors.CalcExactMolWt(molecule_obj))
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
    


def calculate_generations(all_reactants, all_reactions, rxn_ids, num_generations, generation_count):
    print(f"\t{generation_count}...")
    len_reactants_start = len(all_reactants)
    
    # compute list of all possible reaction combinations from global reactants list
    if debug: print("\t\tFinding all possible reactant pair combinations...")
    reactant_combinations = list(itertools.combinations(all_reactants, 2)) # assume always 2 reactants per reaction
    if debug: print("ALL REACTANTS:", all_reactants)
    if debug: wait = input('Press enter to continue or Ctrl+C to exit...')
    if debug: print("ALL REACTANT COMBINATIONS:", reactant_combinations)
    if debug: wait = input('Press enter to continue or Ctrl+C to exit...')
    
    # loop through reactant combinations to call react() and get product data
    if debug: print("\t\tComputing any new reactions and logging information on any new products formed...")
    num_new_prods = 0
    for rxn in reactant_combinations:
        reactant_a = rxn[0]
        reactant_b = rxn[1]
        for rxn_obj in all_reactions:
            rxn_id = rxn_ids[rxn_obj]
            if not reaction_already_computed(reactant_a = reactant_a,
                                             reactant_b = reactant_b,
                                             rxn_id = rxn_id):
                products = react(reactant_a = reactant_a,
                                 reactant_b = reactant_b,
                                 rxn_obj = rxn_obj,
                                 rxn_id = rxn_id)
                if debug: print(products)
                if products != None:
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
                            num_new_prods += 1 # for break condition
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
        # if generation passes user input limit num_generations, break
        return print("All done! Reached maximum number of generations to compute.")
    elif num_new_prods == 0:
        # if no new products are formed this generation, trigger end of script
        return print("All done! Ran out of reactions to do.")
    else:
        calculate_generations(all_reactants = all_reactants,
                              all_reactions = all_reactions,
                              rxn_ids = rxn_ids,
                              num_generations = num_generations, # stop condition never changes from main(), just keep passing to next generation
                              generation_count = generation_count)



def export_network_to_csv_with_self_joins():
    """
    Export "molecule" table in DB to network.csv with timestamp in file name.
    
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


def get_molecule_table():
    df = pd.read_sql("""SELECT * FROM molecule;""",con=engine)
    return df


def export_network_to_csv():
    """
    Export "molecule" table in DB to network.csv with timestamp in file name.
    """
    df = get_molecule_table()
    df.to_csv("output/network_" + timestamp() + ".csv", index=False)



def load_rxn_network_from_db():
    """
    Take Molecule table from db / csv export and load into graph management
    system of choice (i.e. NetworkX or Neo4j).
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
    all_mols_file = open('input/input_molecules.txt','w')
    for mol in all_mols:
        mol = mol.split(r"'")[1]
        # print(mol)
        all_mols_file.write(mol + "\n")
    all_mols_file.close()


def import_molecules():
    # clear molecule table
    
    # load with input/input_molecules.txt
    molecules = codecs.open("input/input_molecules.txt",'r','utf-8').read().split("\n")
    molecules = [mol for mol in molecules if (mol != None) and (mol != "")]
    
    # get molecule objects and create records in db with generation_formed = 0 being parent generation
    mol_objs = [Chem.MolFromSmiles(mol) for mol in molecules]
    for mol in mol_objs:
        create_molecule_record(molecule_obj = mol,
                               generation_count = 0,
                               reactant_1 = None,
                               reactant_2 = None)
    
    # return RDKit molecule objects
    return mol_objs

def import_reactions():
    # clear reaction table
    
    # load with input_reactions.txt
    rxns = codecs.open("input/input_reactions.txt",'r','utf-8').read().split("\n")
    rxns = [rxn for rxn in rxns if (rxn != "") and (rxn != None) and (len(rxn) > 2)]
    
    # create dict of rxn objects and their rxn_id while updating reaction table in db
    all_rxns = []
    rxn_ids = {}
    for smarts_str in rxns:
        print(smarts_str)
        rxn_obj = AllChem.ReactionFromSmarts(smarts_str)
        rxn_ids[rxn_obj] = get_or_create_reaction_record(smarts_str) # create record in db
        all_rxns.append(rxn_obj) # convert to reaction object
    
    # return RDKit reaction objects and rxn_ids
    return (all_rxns, rxn_ids)

def apply_smiles_to_mol(row):
    return Chem.MolFromSmiles(row['smiles_str'])

def apply_exact_mass(row):
    return Chem.rdMolDescriptors.CalcExactMolWt(row['mol_obj'])

def exact_mass_bar_plot():
    # load molecule table
    df = get_molecule_table()
    
    # removed generation formed for input molecules (no reactant FKs available)
    
    
    # get exact mass by molecule
    df['mol_obj'] = df.apply(apply_smiles_to_mol, axis=1)
    df['exact_mass'] = df.apply(apply_exact_mass, axis=1)
    
    # do summary by generation & plot total exact_mass by generation
    df = df.groupby("generation_formed")['exact_mass'].sum()
    print(df)
    df.plot(x='generation_formed', y='exact_mass', kind="bar")
    plt.show()
    


def main():
    """
    Configure and run it!
    """
    # clear all records from db table (generated from last time main() was run)
    print("Preparing database...")
    clear_db_table()
    
    print("Loading input molecules...")
    # define input molecules (SMILES format)
    # test molecules taken from here: http://www.rdkit.org/docs/RDKit_Book.html#reaction-smarts
    # give input molecules generation_formed = 0
    input_molecules = import_molecules()
    # print(input_molecules[0]) # should return __repr__ of molecule object; if loading molecule failed, returns None
    
    print("Loading input reaction rules...")
    # define reactions in SMARTS format
    (input_reactions, rxn_ids) = import_reactions()
    
    
    # compute all generations, up to num_generations (or until products list is empty,
    # i.e. no new reactions.)
    # initialize generation_count as 0; variable will pass through recursively
    # to increment by 1 each generation
    print("Calculating generation number:")
    try:
        calculate_generations(all_reactants = input_molecules,
                              all_reactions = input_reactions,
                              rxn_ids = rxn_ids,
                              num_generations = 5,
                              generation_count = 1)
        export_network_to_csv()
        exact_mass_bar_plot()
    except KeyboardInterrupt:
        export_network_to_csv()
        print("Keyboard Interrupt received. Output exported.")
    # export_network_to_csv_with_self_join() # very slow; just load into network
    
    


# call the function!
main()
# scrape_molecules()
# exact_mass_bar_plot()





