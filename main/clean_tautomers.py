from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

import ambit_tautomer

# molVs' smiles might not be identical to rdkit's so I used RDKit instead of molVs
enum = rdMolStandardize.TautomerEnumerator()

def clean_taut(dg, dg_execute):
	subset = dg_execute.subset
	universe = dg_execute.universe
	# The following line may seem redundant but it is needed for retaining the indices of subset
	graphs = [g for g in subset]
	# list  of smiles associated with each Graph object
	mod_subset_smiles = [g.smiles for g in subset]

	cdk_molecules = [ambit_tautomer.smilesToMolecule(smiles) for smiles in mod_subset_smiles]
	# convert to cdk unique smiles TODO: these smiles are perhaps not canonical.
	cdk_smiles = [ambit_tautomer.smiles_from_molecule(mol) for mol in cdk_molecules]
	# a mapping of the Graph objects with their corresponding DGVertex objects
	mod_dgverts = {v.graph: v for v in dg.vertices if v.graph in subset}
	# a mapping {smiles: tuple(smiles)} between the smiles of a given molecule
	# and those of all possible tautomers for that molecule
	taut_tautset_dict = {}
	for smiles in cdk_smiles:
		taut_tautset_dict[smiles] = ambit_tautomer.generateTautomers(smiles, "CMI")
	tautomer_sets = tuple(taut_tautset_dict.values())
	class_ids = {}
	# create initial class ids
	for i in range(len(tautomer_sets)):
		class_ids[tautomer_sets[i]] = i
	# find sets with common elements and assign them the same class id
	for i in range(len(tautomer_sets)):
		for j in range(i+1, len(tautomer_sets)):
			for item in tautomer_sets[i]:
				if item in tautomer_sets[j]:
					#print(f'Classes {i} and {j} have common elements')
					if class_ids[tautomer_sets[i]] < class_ids[tautomer_sets[j]]:
						class_ids[tautomer_sets[j]] = class_ids[tautomer_sets[i]]
					else:
						class_ids[tautomer_sets[i]] = class_ids[tautomer_sets[j]]
					break

	reversed_dict = {} # {class id: [all relevant smiles in MOD]}
	
	for i in range(len(cdk_smiles)):
		for tautomer_class, class_id in class_ids.items():
			if cdk_smiles[i] in tautomer_class:
				if class_id in reversed_dict.keys():
					reversed_dict[class_id].append(mod_subset_smiles[i])
				elif class_ids[tautomer_class] not in reversed_dict.keys():
					reversed_dict[class_id] = [mod_subset_smiles[i]]
	
	#print("There are {0} unique tautomer classes among {1} molecules".format(len(reversed_dict.keys()), len(subset)))
	
	# a mapping {smiles: smiles}, the one that nees to be kept and the ones that need to be removed
	to_remove = {}

	for seen_tautomers in reversed_dict.values():
		if len(seen_tautomers) > 1:
			for i in range(1, len(seen_tautomers)):
				to_remove[seen_tautomers[0]] = seen_tautomers[i]

	p = GraphPrinter()
	p.simpleCarbons = True
	p.withColour = True
	p.collapseHydrogens = True

	# for each list of redundant tautomer smiles and the corresponding one tautomer to be kept
	for to_keep, item_to_remove in to_remove.items():
		# for each item in the list
		index = mod_subset_smiles.index(item_to_remove)
		index_to_keep = mod_subset_smiles.index(to_keep)
		postSection(f"Removed")
		graphs[index].print(p)
		postSection(f"Kept")
		graphs[index_to_keep].print(p)

		# The DGVertex associated with the molecule to be removed
		dg_vertex = mod_dgverts[graphs[index]]
		# add a fake edge to the preserved graph to account for the removed graph (to avoid loss of reaction info)
		for e in dg_vertex.inEdges:
			for source in e.sources:
				d = Derivations()
				d.left = [source.graph]
				d.rules = e.rules
				d.right = [graphs[index_to_keep]]
				b.addDerivation(d)
				#print(f"Addded fake edge {d}")
		subset.remove(graphs[index])
		universe.remove(graphs[index])
		print(f"Removing {graphs[index]} and keeping {graphs[index_to_keep]}")
	return subset, universe

def clean_taut_rdkit(dg, dg_execute):
	"""
	Clean tautomers from the derivation graph

	Keyword arguments:
	dg          --- a DG (DerivationGraph) object instance holding the network
	dg_execute  --- a DGBuildExecute object instance

	return: the cleaned subset, universe 
	"""
	subset = dg_execute.subset
	universe = dg_execute.universe
	# The following line may seem redundant but it is needed for retaining the indices of subset
	graphs = [g for g in subset]
	# list  of smiles associated with each Graph objectts/reac-space-exp/rules/michaelAddition.py')
	mod_subset_smiles = [g.smiles for g in subset]
	# a mapping of the Graph objects with their corresponding DGVertex objects
	mod_dgverts = {v.graph: v for v in dg.vertices if v.graph in subset}
	# the list of RDKit canonical SMILES
	rdkit_subset_smiles = []
	# a mapping {smiles: tuple(smiles)} between the smiles of a given molecule
	# and those of all possible tautomers for that molecule
	taut_tautset_dict = {}
	# a mapping {smiles: smiles} of the tautomer that needs removing to the one that should be kept
	to_remove = {}
	# Populate all possible tautomers and add to the above empty dict
	for mod_smiles in mod_subset_smiles:
		mol = Chem.MolFromSmiles(mod_smiles)
		# convert into rdkit canonical smiles
		rdkit_smiles = Chem.MolToSmiles(mol)
		rdkit_subset_smiles.append(rdkit_smiles)
		# generate Mol objects of all possible tautomers
		all_tauts = enum.Enumerate(mol)
		# convert into smilesif item in cdk_smiles:
				# TODO: complete
		all_taut_smiles = tuple(Chem.MolToSmiles(taut) for taut in all_tauts)
		taut_tautset_dict[mod_smiles] = all_taut_smiles
	# a tuple of tuples, containing tautomer classes as elements
	list_tautset = tuple(taut_tautset_dict.values())

	# I thought assigning class ids was a decent way of seeing which tautomer classes are identical
	# {tautomer_class: class_id}
	class_ids = {}
	# create initial class ids
	for i in range(len(list_tautset)):
		class_ids.update({list_tautset[i]: i})

	# for debugging purpose
	#ids_tonote = []
	#f = open('incomplete_match.txt', 'w')
	# search for matching tautomer classes
	for i in range(len(list_tautset)):
		for j in range(i+1, len(list_tautset)):
			for item in list_tautset[i]:
				if item in list_tautset[j]:
					#if list_tautset[j] != list_tautset[i]:
						#ids_tonote.append(class_ids[list_tautset[i]])
						#f.write('Classes {0} and {1} have common elements but are not identical'.format(list_tautset[i], list_tautset[j]))
						#f.write('\n')
					#print(f'Classes {i} and {j} have common elements')
					#class_ids[list_tautset[i]] = class_ids[list_tautset[i]]
					class_ids[list_tautset[j]] = class_ids[list_tautset[i]]
					#print('Updated class ids to ', class_ids[list_tautset[i]])
					break
	#f.close()
	# Now, since we're done finding identical class ids
	# A dictionary of observed tautomer classes mapped as {class id: [list of smiles]}
	dict_observed_tauts = {}
	for taut_class, class_id in class_ids.items():
		for i in range(len(rdkit_subset_smiles)):
			if rdkit_subset_smiles[i] in taut_class:
				# note that these are the smiles of the graph generated by mod (not rdkit's smiles)
				if class_id in dict_observed_tauts.keys():
					dict_observed_tauts[class_id].append(mod_subset_smiles[i])
				else :
					dict_observed_tauts[class_id] = [mod_subset_smiles[i]]
	print('{0} unique tautomer classes in {1} molecules'.format(len(dict_observed_tauts.keys())
								, len(subset)))
	# debugging which molecules are causing trouble
	#print(class_ids)
	# add to the dictionary of items to be removed
	for taut_class in dict_observed_tauts.values():
		if len(taut_class) > 1:
			for i in range(1, len(taut_class)):
				to_remove[taut_class[i]] = taut_class[0]
	for smiles in to_remove.keys():
		# The DGVertex associated with the molecule to be removed
		dg_vertex = mod_dgverts[graphs[mod_subset_smiles.index(smiles)]]
		# add a fake edge to the preserved graph to account for the removed graph (to avoid loss of reaction info)
		for e in dg_vertex.inEdges:
			for source in e.sources:
					d = Derivations()
					d.left = [source.graph]
					d.rules = e.rules
					d.right = [graphs[mod_subset_smiles.index(to_remove[smiles])]]
					b.addDerivation(d)
		subset.remove(graphs[mod_subset_smiles.index(smiles)])
		universe.remove(graphs[mod_subset_smiles.index(smiles)])
	return subset, universe