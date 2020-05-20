from rmgpy.rmg.main import RMG
from rdkit import Chem
from rmgpy.species import Species
from rmgpy.solver.simple import SimpleReactor
import networkx as nx
import matplotlib

# This line below is needed when trying to print matplotlib's figures in PyCharm

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

"""
Testing RMG features basically
"""
print('Hello! This is a test.')

# The input has to be in the form of an RMG "job" file
# Refer to the doc for details
file_path = 'input.py'
# an instance of the RMG API
my_rmg = RMG(input_file=file_path)
my_rmg.execute()


def clean_prev(current_file):
    """
    Intended to "clean" (remove) simpleReactor() from an input.py
    :param current_file: the file that was just used for a run
    :return: code without simpleReactor()
    """
    prev_input = open(current_file)
    prev_code = prev_input.readlines()
    prev_input.close()
    # I will be removing simpleReactor() lines that contained in the previous file
    starting_line = 0
    end = 9999
    found = [False, False]
    while not found[1]:
        for line in prev_code:
            if "simpleReactor" in line:
                found[0] = True
                starting_line = prev_code.index(line)
            if found[0] and line == ')':
                found[1] = True
                end = prev_code.index(line)
                break
        if not found[1]:
            raise ValueError('Cannot find simpleReactor(), please check input file')
    del prev_code[starting_line:end + 1]  # the last line needs to be included hence end+1
    return prev_code


# I couldn't find a way to access these from the reaction model itself
# TODO: Access these from the models, instead of providing it explicitly
temperature = 800
pressure = 1.0


def generate_new_input(current_file, new_file, reaction_model, include_edge=False):
    new_input = open(new_file, 'w')
    # remove simpleReactor from prev file
    raw_code = clean_prev(current_file)
    new_code = []
    inerts = ['He', 'Ne', "Ar"]
    for sp in reaction_model.core.species:
        # The inert gases which RMG adds by default must not be saved in the input file
        if sp.label not in inerts:
            new_code.append('\n')
            new_code.append('species(')
            new_code.append("\tlabel='{0}',".format(sp.label))
            new_code.append('\treactive=True,')
            new_code.append(f"\tstructure=SMILES('{sp.smiles}')")
            new_code.append(')')
    # if include_edge:
    # do the same for edge species
    # now add a simpleReactor()

    new_code.append('\n')
    new_code.append('simpleReactor(')
    new_code.append(f"\ttemperature=({temperature},'K'),")
    new_code.append(f"\tpressure=({pressure},'bar'),")
    new_code.append(f"\tterminationTime=(1e0,'s'),")  # by default make all reactions run for 1s
    new_code.append("\tinitialMoleFractions={")

    # I'm assuming we won't add edge species
    # We'll have to expand this method if we do so

    # let all species have equal mole fractions initially
    # # i did a -3 there because core.species contains He, Ne and Ar which we wanna exclude
    average_mol_frac = 1.0 / (len(reaction_model.core.species)-3)
    for sp in reaction_model.core.species:
        if sp.label not in inerts:
            new_code.append(f"\t\t'{sp.label}' : {average_mol_frac},")
            if reaction_model.core.species.index(sp) == (len(reaction_model.core.species) - 1):
                new_code[len(new_code) - 1].replace(',', '')  # remove the comma for the last item
    new_code.append('\t}')
    new_code.append(')')

    # finally, write all of this
    new_input.writelines(raw_code)  # lines from the old input.py don't need an \n added
    for line in new_code:
        new_input.write('\n')
        new_input.write(line)

    new_input.close()


print('Generating input file for next generation')
generate_new_input(file_path, 'input_1.py', my_rmg.reaction_model)
# note that the above process will take long (several minutes)

network_graph = nx.DiGraph()

# create a list of species and reactions for this generation
# Don't forget to unpack the species and reactions lists before adding them into a single list
species_list = [*my_rmg.reaction_model.core.species, *my_rmg.reaction_model.edge.species]
reactions_list = [*my_rmg.reaction_model.core.reactions, *my_rmg.reaction_model.edge.reactions]


def generate_network(species, reactions):
    for sp in species:
        # You can pass the entire Species object too (which has a lot of info)
        # Each species has a unique index integer id during an individual run
        # I passed sp.index for the sake of labels in the graph
        network_graph.add_node(sp.index)
    for reaction in reactions:
        try:
            for reactant, product in reaction.pairs:
                if reaction.is_forward:
                    network_graph.add_edge(reactant.index, product.index)
                else:
                    network_graph.add_edge(product.index, reactant.index)
        except(ValueError):
            # Sometimes it can't find reactant-product pairs, so just let those cases slide
            pass
    # once we're done adding nodes and edges to the graph
    # disable with_labels when passing the entire Species object (which turns into a huge string)
    nx.draw(network_graph, with_labels=True, arrows=True)
    # plt.savefig('glycine-hcn.png')
    plt.show()


generate_network(species_list, reactions_list)

# TODO: find a way to access networks generated by RMG
# as of now it prints an empty list.
print('Printing list of networks')
print(my_rmg.reaction_model.network_list)

## I was trying something here
## The code throws an exception no matter how hard you try
## The methods require RMG#execute() in one way or the other

# read the SDF containing glucose and glycine
# sd_supp = Chem.SDMolSupplier('../data/test_1.sdf')
# for reactant in sd_supp: # Work with only one reactant from the list for now
#    print(Chem.MolToSmiles(reactant))
# reactant_2 = Chem.SDMolSupplier('HCN.sdf')[0] # HCN is the second reactant

# Create an RMG Species object from a smiles string corresponding to
# the reactant
# sp1 = Species(smiles=Chem.MolToSmiles(reactant_1))
# sp2 = Species(smiles=Chem.MolToSmiles(reactant_2))

# my_rmg.load_database()
# performs the reaction
# rxns = reac.react_species((reactant_1, reactant_2))
# print(rxns)
# reactor = SimpleReactor(T=(800,'K'), P = (1e0,'bar'), initial_mole_fractions={'glucose':0.5, 'HCN':0.5})
