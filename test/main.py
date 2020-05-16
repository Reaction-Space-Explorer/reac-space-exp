import rmgpy.rmg.main as mn
from rmgpy.rmg.main import RMG
from rdkit import Chem
from rmgpy.species import Species
from rmgpy.solver.simple import SimpleReactor
import rmgpy.rmg.react as reac

"""
Testing RMG features basically
"""
print('Hello! This is a test.')
# The input has to be in the form of an RMG "job" file
# Refer to the doc for details
file_path = 'input.py'

#an instance of the RMG API
my_rmg = RMG(input_file=file_path)
my_rmg.execute()

# TODO: find a way to make these actually work
# as of now they print empty lists.
print('Printing the list of new species generated {0}'.format(my_rmg.reaction_model.new_species_list))
print('Printing output species list {0}'.format(my_rmg.reaction_model.output_species_list))
print('Printing list of output reactions')
print(my_rmg.reaction_model.output_reaction_list)
print('Printing list of networks')
print(my_rmg.reaction_model.network_list)

## I was trying something here
## The code throws an exception no matter how hard you try
## The methods require RMG#execute() in one way or the other

# read the SDF containing glucose and glycine
#sd_supp = Chem.SDMolSupplier('data/test_1.sdf')
#reactant_1 = sd_supp[0] # Work with only one reactant from the list for now
#print(Chem.MolToSmiles(reactant_1))
#reactant_2 = Chem.SDMolSupplier('HCN.sdf')[0] # HCN is the second reactant

# Create an RMG Species object from a smiles string corresponding to
# the reactant
#sp1 = Species(smiles=Chem.MolToSmiles(reactant_1))
#sp2 = Species(smiles=Chem.MolToSmiles(reactant_2))

#my_rmg.load_database()
# performs the reaction
#rxns = reac.react_species((reactant_1, reactant_2))
#print(rxns)
#reactor = SimpleReactor(T=(800,'K'), P = (1e0,'bar'), initial_mole_fractions={'glucose':0.5, 'HCN':0.5})

#glucose = Species(smiles="C(C1C(C(C(C(O1)O)O)O)O)O")
#hcn = Species(smiles='C#N')

#my_rmg.execute()