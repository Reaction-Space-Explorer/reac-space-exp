import networkx as nx
import itertools
import random
import matplotlib.pyplot as plt


def react(molecules):
    """
    Randomly creating a product from reactants.
    
    This function would be replaced by actual cheminformatics reaction node,
    which receives a list of reactants and outputs a list of products.

    Before sending to this function, send through "can_react()" filter to test
    whether the reactant combination matches all reaction criteria.
    """
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    return random.choice(alphabet)


# initialize graph and global vars
rxn_network = nx.DiGraph()
all_reactants = []
all_rxns = []

def calculate_generation(input_molecules, num_generations, generation_count):
    print(f"Calculating generation {generation_count}...")
    # initialize products for this generation
    # will be added to reactants list for next generation
    all_products = []
    
    # load all_reactants global var with input molecules,
    # as well as loading nodes into graph if not already existing
    for m in input_molecules:
        if m not in all_reactants:
            all_reactants.append(m)
        if not rxn_network.has_node(m):
            rxn_network.add_node(m) # insert kwargs here for species data

    # compute all possible reaction combinations from global reactants list
    reactant_combinations = list(itertools.combinations(input_molecules, 2)) # assume always 2 reactants per reaction
    
    # loop through reactant combinations to call react() and get product data
    rxn_outcomes = []
    for rxn in reactant_combinations:
        if rxn not in all_rxns:
            # if reaction has not been computed and outcome stored, then do reaction
            reactant_a = rxn[0]
            reactant_b = rxn[1]
            product = react([reactant_a, reactant_b])
            if product not in all_products:
                if product not in all_reactants:
                    all_reactants.append(product) # insert product as reactant for next generation
                if not rxn_network.has_node(product):
                    rxn_network.add_node(product) # insert kwargs here for species data
            rxn_data = [reactant_a, reactant_b, product] # insert other data of interest from a reaction here in a data dict {}
            rxn_outcomes.append(rxn_data) # append data on newly processed reactions
            all_rxns.append(rxn_data[:2]) # append meta data only onto all_rxns list (right now no reaction data associated, so only reaction + product species names available)
    
    # push new reaction outcomes to new edges in graph (referencing nodes loaded from above, including previous generations)
    for rxn in rxn_outcomes:
        # add kwargs to edges here for reaction information, i.e. rates,
        # also reaction ID (because multiple reactions can form the same product,
        # and reactants can be used in multiple reactions)
        rxn_network.add_edge(rxn[0],rxn[2]) # add edge between reactant a (ix 0) and product (ix 2)
        rxn_network.add_edge(rxn[1],rxn[2]) # add edge between reactant b (ix 1) and product (ix 2)

    
    # if break condition reached, stop, otherwise run next generation
    generation_count += 1
    if generation_count >= num_generations:
        return print("All done! Reached maximum number of generations to compute.")
    elif not rxn_outcomes:
        return print("All done! Ran out of reactions to do.")
    else:
        calculate_generation(all_reactants, num_generations, generation_count)




# initialize first generation with starter molecules
input_molecules = ['A','B','C','D','E'] # replace with actual input molecule set
num_generations = 5 # initialize maximum number of generations to run

# compute first generation and all generations after, up to num_generations
calculate_generation(input_molecules, num_generations, generation_count = 0)

# visualize reaction network
nx.draw(rxn_network, with_labels=True, arrows=True)
plt.show()








