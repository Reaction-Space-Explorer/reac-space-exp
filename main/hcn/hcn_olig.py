with_formaldehyde = False
include("../main.py")

postChapter("HCN Oligmerisation")

hcn = smiles("C#N", name="Hydrogen Cyanide")
ammonia = smiles("N", name="Ammonia")
water = smiles("O", name="Water")

'''dg = DG.load(inputGraphs, inputRules, "5_dec17.dg")
print("Finished loading from dump file")'''

# Number of generations we want to perform
generations = 6

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
		res = b.execute(addSubset(subset) >> addUniverse(universe) >> strat, verbosity=2)
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
		export_to_neo4j(dg_obj = dg, generation_num = gen)
		write_gen_output(subset, gen+1, reaction_name="hcn_olig")
	print('Completed')


postChapter("Molecules with possible valency issues.")
for v in dg.vertices:
	if "[" in v.graph.smiles:
		print(f"Warning! Possible valency issue in molecule {v.graph.smiles}")
		v.graph.print()
		#find_substruct_producer(dg, v.graph, print_rule=True, condense_str=False)
		#trace_pathways(dg, v.graph.smiles, seed_mol=hcn, exact_molecule=True)


# Dump the dg so it can be loaded again quickly without having to generate it from scratch.
f = dg.dump()
print("Dump file: ", f)

count_rules(dg)

'''with open("smiles5.txt", "w") as wr:
	for v in dg.vertices:
		wr.write(f"{v.graph.smiles}\n")'''
		#if exoaminering.monomorphism(v.graph, labelSettings=
		#	LabelSettings(LabelType.Term, LabelRelation.Specialisation)) > 0:
		#	v.graph.print()
		#	print("Found relevant structure!")

check_sdf_matches(dg, "../../data/HCNfixed.sdf")

iminol = smiles("[N]=[C]O", name="iminol", add=False)
aminol = smiles("[N][C]O", name="aminol", add=False)
diamine = smiles("N[C]N", name="diamine", add=False)
find_substruct_producer(dg, iminol, print_rule=True)
find_substruct_producer(dg, aminol, print_rule=True)
#find_substruct_producer(dg, imine, print_rule=True)
##find_substruct_producer(dg, enol, print_rule=True)
#find_substruct_producer(dg, diol, print_rule=True)
find_substruct_producer(dg, diamine, print_rule=True)

rule_names= ["nitriles", "Ring Closure", "Exoamine", "Ammonolysis", "CN", "Amine", "Amide"]
# see what rules containing these words are doing
#for name in rule_names:
#	print_reaction(dg, reac_name=name)