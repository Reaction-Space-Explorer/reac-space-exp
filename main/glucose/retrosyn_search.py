"""

There's a fraction of test set species that couldn't be matched within the first 5 gens and we couldn't go beyond that
due to computational resource limitations. We try to run "retrosynthetic" networks using the same set of reaction rules
(most of which support a backwards reaction, making retrosynthesis possible) to anticipate how long formation of these
species would take.

"""

import time

#from rdkit.Chem import MolFromSmiles, MolToSmiles, Kekulize
# Formaldehyde was produced in G1 or G2 of glucose degradation itself, so it should be available.
# Thus we can use the formaldehyde rule for Cannizarro

with_formaldehyde = True

include("../main.py")
postChapter('Retrosynthesis to find additional test set matches')

# The following includes whatever is in the test set SDF that couldn't be matched by G5
structures_to_search = ['CCC(=O)O',
				'OC1=CC=CC=C1',
				'OC1=C(O)C=CC=C1',
				'CC1CCC(=O)C1=O',
				'CC1=CC=CC(O)=C1O',
				'CC1=CC(O)=C(O)C=C1',
				'CC1=C(O)C=CC(O)=C1',
				'CC1CC(=O)C(C)C1=O',
				'CC1CC(C)C(=O)C1=O',
				'O=CC1=C(O)C=C(O)C=C1',
				'COC1=C(C)C(O)=CC=C1',
				'CC(=O)C1=C(O)C(O)=CC=C1',
				'CC(=O)C1=C(O)C=CC(O)=C1']

# starting reactants.
form = smiles("C=O", name="Formaldehyde")
water = smiles("O", name="Water")

# TODO: find out whether the below rule is actually needed for the purpose or what I used it for.
# For testing something
ether_cleavage = ruleGMLString("""rule [
	ruleID "Cleave ether (test rule)"
	left [
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 3 target 4 label "-" ]
	]
	context [
		node [ id 1 label "H" ]
		node [ id 2 label "C" ]
		node [ id 3 label "C" ]
		node [ id 4 label "O" ]
		edge [ source 4 target 5 label "-" ]
		node [ id 5 label "C" ]
	]
	right [
		edge [ source 1 target 4 label "-" ]
		edge [ source 2 target 3 label "=" ]
	]
]
""")
ether_cleavage.print()

# Read the output of the network
smiles_gen_map = {}
smiles_upto_g5 = []
with open("glucose_degradation_output_10mar.txt") as gen_out:
	for line in gen_out.readlines():
		comps = line.split("\t")
		gen = comps[0][1]
		# Get rid of new-line character
		smi = comps[1].replace("\n", "")
		smiles_gen_map[smi] = gen
		smiles_upto_g5.append(smi)
	print("Done reading output")

'''for s in all_smiles:
	temp = smiles(s, add=False)
	for fb in forbidden:
		if fb.monomorphism(temp, labelSettings=
			LabelSettings(LabelType.Term, LabelRelation.Specialisation)) != 0:
			print(f"Printing {temp}")
			temp.print(p)'''

# Number of generations we want to perform
generations = 4
#find_substruct_producer(dg, ethoxy_ether)

#propanoic = smiles("CCC(=O)O", name="Propanoic acid", add=False)
#ribitol = smiles("OC(O)C(O)C(O)C(O)CO", name="Ribitol")

ignore = ["C(O)=O", "CO"] # stuff produced by formaldehyde itself
def check_for_precursor(subset, retro_gen):
	product_smiles_list = [g.smiles for g in subset]
	for item, gen in smiles_gen_map.items():
		# if any molecule in the network upto G5 shows up in this product suite, then we know a precusor
		# exists in the forward running glucose network
		if item in product_smiles_list and item not in ignore:
			print(f"Found {item} as a possible precusor! Took {retro_gen} generations of retrosynthesis")
			subset[product_smiles_list.index(item)].print(p)


# Generate a full fledged network starting with each of the test-set molecules we're still trying to find
for smi in structures_to_search:
	postSection(f"Searching for {smi}")
	print('Searching for {smi}')
	to_search = smiles(smi, name="Test set item")
	to_search.print(p)
	dg = DG(graphDatabase=inputGraphs,
		labelSettings=LabelSettings(LabelType.Term, LabelRelation.Specialisation))
	subset = inputGraphs
	universe = []
	postSection("Possible Precusors")
	with dg.build() as b:
		for gen in range(generations):
			postSection(f"Generation {gen+1}")
			start_time = time.time()
			print(f"Starting round {gen+1}")
			res = b.execute(addSubset(subset) >> addUniverse(universe) >> strat, verbosity=2)
			end_time = time.time()
			print(f"Took {end_time - start_time} seconds to complete round {gen+1}")
			print('Products subset size:', len(res.subset))
			subset, universe = res.subset, res.universe
			check_for_precursor(subset, gen+1)
	d = dg.dump()
	print("Dump file:", d)
	# now carefully remove this starting reactant before adding a searching again for a diff one
	inputGraphs.remove(to_search)
			
'''	for g in subset:
				if g.isomorphism(propanoic) == 1:
					print("found propanoic acid!")
		print('Completed')'''
