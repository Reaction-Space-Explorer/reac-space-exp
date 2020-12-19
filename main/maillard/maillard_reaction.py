with_formaldehyde = False

include("../main.py")

postChapter("Maillard Reaction")

glycine = smiles("NCC(=O)O", name="Glycine")
open_glucose = smiles("O=CC(O)C(O)C(O)C(O)C(O)", "Open Chain Glucose")
water = smiles("O", name="Water")

'''dg = DG.load(inputGraphs, inputRules, "3_dec15.dg")
print("Finished loading from dump file")'''

# Number of generations we want to perform
generations = 4

dg = DG(graphDatabase=inputGraphs,
	labelSettings=LabelSettings(LabelType.Term, LabelRelation.Specialisation))

subset = inputGraphs
universe = []
# In the following block, apart from generating the reactions, we may print structures
# and reactions forming them that are not in the MS
#postSection("Structures not found in MS")
with dg.build() as b:
	for gen in range(generations):
		start_time = time.time()
		print(f"Starting round {gen+1}")
		res = b.execute(addSubset(subset) >> addUniverse(universe) >> strat, verbosity=8)
		end_time = time.time()
		print(f"Took {end_time - start_time} seconds to complete round {gen+1}")
		print(f'Products in round {gen+1}: {len(res.subset)}')

		# The returned subset and universe do not contain redundant tautomers
		#subset, universe = clean_taut(dg, res, algorithm="CMI")
		subset, universe = res.subset, res.universe
		#print('Subset size after removal:', len(subset))
		# This step replaces the previous subset (containing tautomers) with the cleaned subset
		#res = b.execute(addSubset(subset) >> addUniverse(universe))
		# now compare how many of these simulations were found in the MS data.
		#export_to_neo4j(dg_obj = dg, generation_num = gen)
		write_gen_output(subset, gen+1, reaction_name="maillard")
	print('Completed')


aminol = smiles("O[C][N]", name="Aminol", add=False)
imine = smiles("[C]=[N]", name="Imine", add=False)
enol = smiles("[C]=[C]O", name="enol", add=False)
diol = smiles("O[C]O", name="diol", add=False)
diamine = smiles("[N][C][N]", name="diamine", add=False)
find_substruct_producer(dg, aminol, print_rule=True)
find_substruct_producer(dg, imine, print_rule=True)
find_substruct_producer(dg, enol, print_rule=True)
find_substruct_producer(dg, diol, print_rule=True)
find_substruct_producer(dg, diamine, print_rule=True)

#print_reaction(dg, "Deamination", print_max=300)
#print_reaction(dg, "of Amines", print_max=300)
#print_reaction(dg, "Enamine Hydration", print_max=300)

# Dump the dg so it can be loaded again quickly without having to generate it from scratch.
f = dg.dump()
print("Dump file: ", f)

with open("maillard_3r_dec15.txt", "w") as wr:
	for v in dg.vertices:
		wr.write(f"{v.graph.smiles}\n")

check_sdf_matches(dg, "../../data/Maillard2.sdf")

count_rules(dg)