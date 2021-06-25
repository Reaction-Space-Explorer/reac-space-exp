with_formaldehyde = False

include("../main.py")
postChapter("HCN Oligmerisation")

hcn = smiles("C#N", name="Hydrogen Cyanide")
ammonia = smiles("N", name="Ammonia")
water = smiles("O", name="Water")

'''dg = DG.load(inputGraphs, inputRules, "4rounds_dec25.dg")
print("Finished loading from dump file")'''

# Number of generations we want to perform
generations = 7

dg = DG(graphDatabase=inputGraphs,
	labelSettings=LabelSettings(LabelType.Term, LabelRelation.Specialisation))

subset = inputGraphs
universe = []
# dump initial reactants as part of "G0"
write_gen_output(subset, generation=0, reaction_name="hcn_olig")

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
		# export some info for further analysis
		export_to_neo4j(dg_obj = dg, generation_num = gen+1)
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

# NOTE: turn this on when you need to export it.
#count_rules_by_gen(dg, 'hcn_olig_output.txt')
#dg.print()


check_sdf_matches(dg, "../../data/HCNfixed.sdf")

# Check if any reactions are producing unwanted tautomers

iminol = smiles("[N]=[C]O", name="iminol", add=False)
aminol = smiles("[N][C]O", name="aminol", add=False)
diamine = smiles("N[C]N", name="diamine", add=False)
find_substruct_producer(dg, iminol, print_rule=True)
find_substruct_producer(dg, aminol, print_rule=True)
#find_substruct_producer(dg, imine, print_rule=True)
##find_substruct_producer(dg, enol, print_rule=True)
#find_substruct_producer(dg, diol, print_rule=True)
find_substruct_producer(dg, diamine, print_rule=True)