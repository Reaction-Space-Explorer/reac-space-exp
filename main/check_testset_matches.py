import pandas as pd

"""

There's a method for checking test set matches in main.py which can do this at the end of reaction
generation itself, but that doesn't take generation into account and
I wanted to see which test set molecule matched in which generation from the output
dump files.

"""
from rdkit.Chem import SDMolSupplier, MolFromSmiles, MolToSmiles, Kekulize

test_set = [mol for mol in SDMolSupplier('../data/NewAlkalineHydrolysisStructures.sdf')]
# Kekulize each species in the test set
for mol in test_set:
	Kekulize(mol)

matched_flag_list = [False for mol in test_set]

match_gen_map = {}

def check_match(candidate_smiles, gen):
	"""
	If A is a substructure of B and B is a substructure of A then the two must be isomorphic. Does that make sense?
	That's the logic I used!
	"""
	candidate = MolFromSmiles(candidate_smiles)
	for m in test_set:
		if m.HasSubstructMatch(candidate) and candidate.HasSubstructMatch(m):
			print(f'{candidate_smiles} matched in {gen}')
			match_gen_map[candidate_smiles] = int(gen[1])
			# update flag for the matched str.
			matched_flag_list[test_set.index(m)] = True


# Open the .txt output
output_data = pd.read_csv('glucose/glucose_degradation_output_10mar.txt',
				 sep='\t', names=['Generation', 'SMILES'])

for i in range(len(output_data)):
	check_match(output_data['SMILES'].iloc[i], output_data['Generation'].iloc[i])


# list molecules that couldn't be matched
print('The following test set molecules didn\'t have a match')
for i in range(len(test_set)):
	if matched_flag_list[i] == False:
		smi = MolToSmiles(test_set[i], kekuleSmiles=True, isomericSmiles=False) #isomericSmiles includes stereochem info
		print(smi)