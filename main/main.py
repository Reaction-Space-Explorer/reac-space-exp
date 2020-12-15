import os
import time
from rdkit.Chem import SDMolSupplier, MolFromSmiles, MolToSmiles, Kekulize

include(os.path.join('..', 'rules/all.py'))
include("mod_to_neo4j_exporter.py")


p = GraphPrinter()
p.simpleCarbons = True
p.withColour = True
p.collapseHydrogens = True

# Forbidden substructures (don't add them to inputGraphs)
# Three and four membeered rings are unstable, any atom with two double bonds is forbidden
# Note: the latter means that even O=C=O (carbon dioxide) will also be forbidden

three_memb = graphGMLString("""graph [
		node [ id 0 label "*" ]
		node [ id 1 label "*" ]
		node [ id 2 label "*" ]
		edge [ source 0 target 1 label "*" ]
		edge [ source 1 target 2 label "*" ]
		edge [ source 2 target 0 label "*" ]
	]""", name= "Three membered ring", add=False)

four_memb = graphGMLString("""graph [
	node [ id 0 label "*" ]
	node [ id 1 label "*" ]
	node [ id 2 label "*" ]
	node [ id 3 label "*" ]
	edge [ source 0 target 1 label "*" ]
	edge [ source 1 target 2 label "*" ]
	edge [ source 2 target 3 label "*" ]
	edge [ source 3 target 0 label "*" ]
]""", name="Four membered ring", add=False)

# some substructures had high electron density
bad_bicycle_1 = smiles("[C]1[O][C](O)[O][C]1", add=False)
bad_bicycle_2 = smiles("[C]1(O)[O][C][O][C]1", add=False)
forbidden = [three_memb, four_memb, smiles('[*]=[*]=[*]', name="Two double bonds", add=False),
			bad_bicycle_1, bad_bicycle_2]

# Load the library of bad ring and aromatic structures.
bad_rings_aromatics = ['BadRingsList.txt', 'BadAromaticsList.txt']
for file_name in bad_rings_aromatics:
	with open(f"../../data/{file_name}") as rings_list:
		list_rings = rings_list.readlines()
		for item in list_rings:
			ring_smiles = item.replace("\n", "")
			ring = smiles(ring_smiles, add=False)
			forbidden.append(ring)


# A list of things that might be forbidden by the library above but is stable
# and should be produced
allowed_structs = [smiles("O=C=O", name="Carbon Dioxide", add=False),
					smiles("[N]=C=O", name="isocyanates", add=False)]

# Make this a global variable (advantage being it can be changed in other reactions individually)
max_mass_limit = 200

def pred(derivation):
	"""
	Keyword arguments:
	d --- a derivation graph object
	"""
	for g in derivation.right:
		for allowed in allowed_structs:
			if g.monomorphism(allowed) > 0:
				return True
		# Allow masses only lower than a certain maximum
		if g.exactMass >= max_mass_limit:
			return False
		for fb in forbidden:
			if fb.monomorphism(g, labelSettings=
			LabelSettings(LabelType.Term, LabelRelation.Specialisation)) > 0:
				return False
		#print(g)
	return True

strat = (
	addSubset(inputGraphs)
	>> rightPredicate[
		pred
		#lambda derivation: all (g.exactMass <= 500
		# and (fb.monomorphism(g) == 0 for fb in forbidden) for g in derivation.right)
	] (inputRules)#repeat[1]()
)

# Use the following when applying a strategy to a DG loaded from a dump file
'''start_time = 0
with dg_new.build() as b:
	for v in dg.vertices:
		if v.graph not in inputGraphs:
			#print(f"Adding {v.graph}")
			inputGraphs.append(v.graph)
		else:
			pass
	# add everything in dg to dg_new
	for e in dg.edges:
		sources = [source.graph for source in e.sources]
		targets = [target.graph for target in e.targets]
		d = Derivations()
		d.left = sources
		d.rules = e.rules
		d.right = targets
		b.addDerivation(d)
	start_time = time.time()
	res = b.execute(addSubset(inputGraphs) >> rightPredicate[pred](rules_list), verbosity=8)
	print("Number of products this round: ", len(res.subset))
end_time = time.time()
print("Time taken for this 4th round: ", end_time-start_time)'''


def check_sdf_matches(dg, sdf_file, draw_structures=True, print_unmatching=False):
	"""
	After generating the network, try to see if any structures match with those in SDF files
	These files were usually created manually, storing structures reported in experimental
	studies. The purpose is to match our simulations with experiments.

	Keyword arguments:
	dg			-- the derivation graph of the network
	sdf_file	-- path to the SDF file
	draw_structures -- whether or not to print the structures in the summary pdf
	"""
	matching_structs = []
	not_matching = []
	postSection('Matching Structures')
	print(f"Checking for matches with structures in {sdf_file}")
	sdfile = SDMolSupplier(sdf_file)
	for mol in sdfile:
		Kekulize(mol)
		smi = MolToSmiles(mol, kekuleSmiles=True)
		mol_graph = smiles(smi, add=False)
		for v in dg.vertices: #dg_new.vertices
			if v.graph.isomorphism(mol_graph) == 1:
				matching_structs.append(mol_graph)
				print(f"Structure {v.graph.smiles} in the network matches a test set molecule!")
		if mol_graph not in matching_structs:
			not_matching.append(mol_graph)
	if draw_structures == True:
		for g in matching_structs:
			g.print(p)
	if print_unmatching == True:
		postSection("Structures not matched yet")
		for g in not_matching:
			g.print(p)

	print(f"{len(matching_structs)} of {len(sdfile)} ({100* len(matching_structs)/len(sdfile)}%)  total structures in the SDF are in the reaction network.")


def write_gen_output(subset, generation, reaction_name):
	"""
	Create an output.txt while storing SMILES and generation number, so that we can later
	compare generation by generation for common structures in different reactions
	"""
	with open(f"{reaction_name}_output.txt", "a") as f:
		for graph in subset:
			f.write(f"G{generation}\t{graph.smiles}\n")


# Count number of times a rule appears in the edges' attributes.
def count_rules(dg):
	"""
	Counts the number of edges associated with each rule in the network and prints the count.
	Keyword arguments:
	dg -- the derivation graph of the network
	"""
	rules_count = []
	for e in dg.edges:
		for rule in e.rules:
			rules_count.append(rule.name)
	rules_used = dict({rule:True for rule in rules_count})
	for rule in rules_used.keys():
		print(f"{rule} reaction count: {rules_count.count(rule)}")

dgprint = DGPrinter()
#dgprint.withRuleName = True
dgprint.withShortcutEdges = True

rp = GraphPrinter() # Think of it as a "Rule Printer"
rp.withIndex = True
rp.withColour = True
# Track what is producing diols
gem_diol = smiles("O[C]O", name="gem diol substruct", add=False)

# Track what's producing methoxy/ethoxy ethers
methoxy = smiles("[C][O]C", name="methoxy ether substruct", add=False)
ethoxy = smiles("[C][O]CC", name="Ethoxy Ether Substruct", add=False)


def find_substruct_producer(dg, substruct, print_max=50, print_rule=False, condense_str=True
							,exact_struct=False):
	"""
	Find reactions producing a certain substructure.
	"""
	postSection(f"{substruct.name} producing reactions")
	count = 0
	for e in dg.edges:
		if count > print_max:
			break
		else:
			for rule in e.rules:
				sources = [source.graph for source in e.sources]
				targets = [target.graph for target in e.targets]
				for g in targets:
					if substruct.monomorphism(g) > 0:
						count += 1
						print(f"Found {substruct.name}! Made by {rule}")
						dg2 = DG(graphDatabase=inputGraphs)
						if condense_str == False:
							dgprint.graphPrinter.collapseHydrogens = False
							dgprint.graphPrinter.simpleCarbons = False
						with dg2.build() as b:
							d = Derivations()
							d.left = sources
							d.rules = [rule]
							d.right = targets
							fake_edge = b.addDerivation(d)
							print("Printing reaction: ", fake_edge)
						dg2.print(dgprint)
						if print_rule == True:
							rule.print(rp)


def print_reaction(dg, reac_name, print_max=50, print_rule=False):
	"""
	Print individual reactions for testing purposes.
	dg -- the DG of the whole network
	name -- name of the reaction (ruleID)
	max -- the maximum number of such reactions to print.
	"""
	postSection(f"{reac_name} Reactions")
	count = 0
	for e in dg.edges:
		if count > print_max:
			break
		else:
			for rule in e.rules:
				if reac_name in rule.name:
					dg2 = DG(graphDatabase=inputGraphs)
					sources = [source.graph for source in e.sources]
					targets = [target.graph for target in e.targets]
					with dg2.build() as b:
						d = Derivations()
						d.left = sources
						d.rules = [rule]
						d.right = targets
						fake_edge = b.addDerivation(d)
						print("Printing reaction: ", fake_edge)
					dg2.print(dgprint)
					if print_rule == True:
						rule.print(rp)
					count += 1


def trace_pathways(dg, mol_smiles, seed_mol, exact_molecule=False):
	"""
	Method to trace the entire molecular pathway producing a certain species
	beginning from the starting substrates to the target species itself

	TODO: Highlight the target molecule in a diff colour. Colour precursors by generation maybe?

	Keyword arguments:
	dg -- the DG object of the network we're looking at 
	mol_smiles -- SMILES of the target molecule
	seed_mol -- the initial reactant which was used as a "seed" when the reaction network 
			was generated
	exact_molecule -- True if the mol_smiles is the exact structure to search for, otherwise
			if False, mol_smiles will be treated as a substructure and pathways for all molecules
			containing those substructures will be printed
	"""
	# molecule being searched for
	target_mol = smiles(mol_smiles, add=False)
	target_indices = [] # indices of the target molecule(s) in the dg's list of vertices (species)
	all_species = [v for v in dg.vertices]
	all_species_graphs = [v.graph for v in all_species]
	if exact_molecule == False:
		for sp in all_species_graphs:
			if target_mol.monomorphism(sp) >= 1: # if it's present as a substructure
				target_indices.append(all_species_graphs.index(sp))
	else:
		for sp in all_species_graphs:
			if target_mol.isomorphism(sp) == 1: # needs to be an exact match
				target_indices.append(all_species_graphs.index(sp))
	# now, for the molecule at the target index, trace precursors
	for index in target_indices:
		target_dg_vert = all_species[index] # the DGVertex associated with the target molecule
		pathway_dg = DG(graphDatabase=inputGraphs)
		sources = []
		targets = [target_dg_vert]
		trace_complete = False
		with pathway_dg.build() as build_pathway:
			while not trace_complete:
				#print(f"Targets: {[t.graph.smiles for t in targets]}")
				next_targets = []
				for targ in targets:
					if targ.graph.isomorphism(seed_mol) == 1:
						trace_complete = True
						break
					else:
						#print(f"Current target: {targ.graph.smiles}")
						for e in targ.inEdges:
							#print(f"working with edge {e.id}")
							sources = [source.graph for source in e.sources]
							next_targets.extend(source for source in e.sources 
									if source not in next_targets)
							for educt in sources:
								d = Derivations()
								d.left = [educt]
								d.rules = [r for r in e.rules]
								d.right = [targ.graph]
								edge = build_pathway.addDerivation(d)
								#print(f"Adding reaction {d} to pathway")
				#print(f"Targets were {targets}")
				targets = next_targets
				#print(f"Targets are {targets}")
		pathway_dg.print()
