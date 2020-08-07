import os
import time
from rdkit.Chem import SDMolSupplier, MolToSmiles

include(os.path.join('..', 'rules/all.py'))

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

forbidden = [three_memb, four_memb, smiles('[*]=[*]=[*]', name="Two double bonds", add=False)]

def pred(derivation):
	"""
	Keyword arguments:
	d --- a derivation graph object
	"""
	for g in derivation.right:
		# Allow masses only lower than a certain maximum
		if g.exactMass >= 200:
			return False
		for fb in forbidden:
			if fb.monomorphism(g, labelSettings=
			LabelSettings(LabelType.Term, LabelRelation.Specialisation)) > 0:
				return False
		#print(g)
	return True

# store the rules in a separate list
rules_list = []
for rule in inputRules:
	rules_list.append(rule)
# now load Cannizarro2 into inputRules so that loading the DG doesn't give an error (rule not found)
# because the dumped dg had Cannizarro 2 reactions in it.
include(os.path.join("..", "rules/cannizarro2.py"))

strat = (
	addSubset(inputGraphs)
	>> rightPredicate[
		pred
		#lambda derivation: all (g.exactMass <= 500
		# and (fb.monomorphism(g) == 0 for fb in forbidden) for g in derivation.right)
	] (inputRules)#repeat[1]()
)

p = GraphPrinter()
p.simpleCarbons = True
p.withColour = True
p.collapseHydrogens = True

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


def check_sdf_matches(dg, sdf_file, draw_structures=True):
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
		print(mol)
		smi = MolToSmiles(mol)
		#smi.replace("", "")
		mol_graph = smiles(smi, add=False)
		for v in dg.vertices: #dg_new.vertices
			if v.graph.isomorphism(mol_graph) == 1:
				matching_structs.append(mol_graph)
				print("Structure {0} of the SDF found in the network!".format(mol_graph.smiles))
				if draw_structures == True:
					v.graph.print(p)
			else:
				not_matching.append(mol_graph)
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
	rules_count = []
	for e in dg.edges:
		for rule in e.rules:
			rules_count.append(rule.name)
	rules_used = dict({rule:True for rule in rules_count})
	for rule in rules_used.keys():
		print(f"{rule} reaction count: {rules_count.count(rule)}")

# Print reactions of just the one reaction alone, in each in its own DG
#to_print = ['Cannizarro', 'Knoevenagel, C, C(=O)A, OR, H |, HC', 'Ester Formation Hydrolysis Exchange',
#			'Knoevenagel, C, C(=O)A, OR, H |, CC', 'Ring Closure']

dgprint = DGPrinter()
dgprint.withRuleName = True
dgprint.withShortcutEdges = True

# Track what is producing diols
'''diol = graphGMLString("""graph [
	node [ id 0 label "C" ]
	node [ id 1 label "O" ]
	node [ id 2 label "O" ]
	node [ id 3 label "H" ]
	node [ id 4 label "H" ]
	edge [ source 0 target 1 label "-" ]
	edge [ source 0 target 2 label "-" ]
	edge [ source 1 target 3 label "-" ]
	edge [ source 2 target 4 label "-" ]
]""", add=False)'''

# Track what's producing methoxy ethers
methoxy_ether = smiles("[O]C", add=False)


#postSection("Methoxy ether producing reactions")
#for item_to_print in to_print:
#count = 0
#postSection(f"{item_to_print} reactions")
# print all reactions
'''for e in dg.edges:
		# Don't print more than 35 of any category
	if count > 50:
		#pass
		break
	else:
		for rule in e.rules:
			if "" in rule.name:
				dg2 = DG(graphDatabase=inputGraphs)
				sources = [source.graph for source in e.sources]
				targets = [target.graph for target in e.targets]
				for g in targets:
					if methoxy_ether.monomorphism(g) > 0:
						count += 1
						print("Found methoxy ether!")
						with dg2.build() as b:
							d = Derivations()
							d.left = sources
							d.rules = [rule]
							d.right = targets
							fake_edge = b.addDerivation(d)
							print("Printing reaction: ", fake_edge)
							#rule.print()
						dg2.print(dgprint)'''