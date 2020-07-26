import os
import time
from rdkit.Chem import SDMolSupplier, MolToSmiles

include(os.path.abspath(os.path.join('..', 'rules/all.py')))
# including these files using MOD's include() so that MOD's functions are callable in them
include("compare_ms.py")
#include('clean_tautomers.py')

#postChapter('Alkaline Glucose Degradation')
postChapter("Formose Reaction")
# starting molecule
#glucose = smiles('OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O', name='Glucose')
formaldehyde = smiles("C=O", name="Formaldehyde")
glycoladehyde = smiles("OCC=O", name="Glycolaldehyde")
#open_glucose = smiles("O=CC(O)C(O)C(O)C(O)C(O)", "Open Chain Glucose")
water = smiles("O", name="Water")

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

strat = (
	addSubset(inputGraphs)
	>> rightPredicate[
		pred
		#lambda derivation: all (g.exactMass <= 500
		# and (fb.monomorphism(g) == 0 for fb in forbidden) for g in derivation.right)
	] (inputRules)#repeat[1]()
)

# Number of generations we want to perform
generations = 3

#postSection('Final Network')
dg = DG(graphDatabase=inputGraphs,
	labelSettings=LabelSettings(LabelType.Term, LabelRelation.Specialisation))

# store the rules in a separate list
rules_list = []
for rule in inputRules:
	rules_list.append(rule)
# now load Cannizarro2 into inputRules so that loading the DG doesn't give an error (rule not found)
# because the dumped dg had Cannizarro 2 reactions in it.
include(os.path.join("..", "rules/cannizarro2.py"))

p = GraphPrinter()
p.simpleCarbons = True
p.withColour = True
p.collapseHydrogens = True

dg = dgDump(inputGraphs, inputRules, "round4_without_cann.dg")
print("Finished loading from dump file")

postSection("Four Membered Rings")
print("Checking for four membered rings")

dg_new = DG(graphDatabase=inputGraphs,
	labelSettings=LabelSettings(LabelType.Term, LabelRelation.Specialisation))

start_time = 0
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
	# Confirming that everything in dg was added to dg_new
	print("Checking for differences")
	diffDGs(dg, dg_new)
	print("If nothing was printed, the two DGs are identical")
	#once done, now we can do our stuff
	start_time = time.time()
	res = b.execute(addSubset(inputGraphs) >> rightPredicate[pred](rules_list), verbosity=8)
	print("Number of products this round: ", len(res.subset))
end_time = time.time()
print("Time taken for this 4th round: ", end_time-start_time)


subset = inputGraphs
universe = []

# In the following block, apart from generating the reactions, we may print structures
# and reactions forming them that are not in the MS
#postSection("Structures not found in MS")
'''with dg.build() as b:
	for gen in range(generations):
		start_time = time.time()
		print(f"Starting round {gen+1}")
		res = b.execute(addSubset(subset) >> addUniverse(universe) >> strat, verbosity=8)
		end_time = time.time()
		print(f"Took {end_time - start_time} seconds to complete round {gen+1}")
		print('Original subset size:', len(res.subset))

		# The returned subset and universe do not contain redundant tautomers
		#subset, universe = clean_taut(dg, res, algorithm="CMI")
		subset, universe = res.subset, res.universe
		#print('Subset size after removal:', len(subset))
		# This step replaces the previous subset (containing tautomers) with the cleaned subset
		#res = b.execute(addSubset(subset) >> addUniverse(universe))
		# now compare how man
		# y of these simulations were found in the MS data.
		#compare_sims(dg, gen+1, print_extra=False)
	print('Completed')'''

'''with open("formose_reaction_4_rounds.txt", "w") as fish:
	for v in dg_new.vertices:
		fish.write(v.graph.smiles)
		fish.write('\n')'''

# compare structures with what the Y&M paper has
sdfile = SDMolSupplier(os.path.join("..", "data/NewAlkalineHydrolysisStructures.sdf"))

matching_structs = []
not_matching = []
'''postSection('Matching Structures')

print("Checking for matches with Y&M's structures")
for mol in sdfile:
	smi = MolToSmiles(mol)
	#smi.replace("", "")
	mol_graph = smiles(smi, add=False)
	for v in dg_new.vertices: #dg.vertices
		if v.graph.isomorphism(mol_graph) == 1:
			matching_structs.append(mol_graph)
			print("Structure {0} of the SDF found in the network!".format(mol_graph.smiles))
			v.graph.print(p)
		else:
			not_matching.append(mol_graph)

print(f"{len(matching_structs)} of {len(sdfile)} ({100* len(matching_structs)/len(sdfile)}%)  total structures in the SDF are in the reaction network.")
'''
postSection("Molecules with possible incomplete valencies")
#postSection("All vertices")
for v in dg.vertices:
	if '[C' in v.graph.smiles:
		v.graph.print(p)

f = dg.dump()
print("Dump file: ", f)

rules_count = []
for e in dg.edges:
	for rule in e.rules:
		rules_count.append(rule.name)

rules_used = dict({rule:True for rule in rules_count})
for rule in rules_used.keys():
	print(f"{rule} reaction count: {rules_count.count(rule)}")
#print("Aldol condensation reaction count: ", rules_count.count("Aldol Condensation"))

#print("Rules used: {0}".format(dict({rule:True for rule in rules_count}).keys()))
# Make a mass spectra (a histogram of the masses) of the molecules
#compare_ms.make_mass_spectra([v.graph.smiles for v in dg.vertices])

# Compare structures with 

# Print reactions of just the one reaction alone, in separate DGs
to_print = ['Cannizarro', 'Knoevenagel, C, C(=O)A, OR, H |, HC', 'Ester Formation Hydrolysis Exchange',
			'Knoevenagel, C, C(=O)A, OR, H |, CC', 'Ring Closure']

dgprint = DGPrinter()
dgprint.withRuleName = True
dgprint.withShortcutEdges = True

postSection("Reactions")
#for item_to_print in to_print:
count = 0
#postSection(f"{item_to_print} reactions")
# print all reactions
'''for e in dg.edges:
		# Don't print more than 35 of any category
#        if count > 35:
#            pass
		#    break
#        else:
	for rule in e.rules:
		if "Aldol" in rule.name:
#        count += 1
			dg2 = DG(graphDatabase=inputGraphs)
			with dg2.build() as b:
				d = Derivations()
				sources = [source.graph for source in e.sources]
				targets = [target.graph for target in e.targets]
				d.left = sources
				d.rules = [rule]
				d.right = targets
				fake_edge = b.addDerivation(d)
				print("Printing reaction: ", fake_edge)
			dg2.print(dgprint)'''
# dump smiles