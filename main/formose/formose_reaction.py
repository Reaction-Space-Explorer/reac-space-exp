# This flag will be used to toggle one pair of extra Cannizarro 2 rules
with_formaldehyde = True

include("../main.py")

postChapter("Formose Reaction")
formaldehyde = smiles("C=O", name="Formaldehyde")
glycoladehyde = smiles("OCC=O", name="Glycolaldehyde")
water = smiles("O", name="Water")

'''dg = dgDump(inputGraphs, inputRules, "../dumps/formose_4_rounds.dg")
print("Finished loading from dump file")'''

# Number of generations we want to perform
generations = 5

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
		print('Original subset size:', len(res.subset))

		# The returned subset and universe do not contain redundant tautomers
		#subset, universe = clean_taut(dg, res, algorithm="CMI")
		subset, universe = res.subset, res.universe
		#print('Subset size after removal:', len(subset))
		# This step replaces the previous subset (containing tautomers) with the cleaned subset
		#res = b.execute(addSubset(subset) >> addUniverse(universe))
		# now compare how many of these simulations were found in the MS data.
		#compare_sims(dg, gen+1, print_extra=False)
		#export_to_neo4j(dg_obj = dg, generation_num = gen)
		write_gen_output(subset, gen+1, reaction_name="formose")
	print('Completed')

# Dump the dg so it can be loaded again quickly without having to generate it from scratch.
f = dg.dump()
print("Dump file: ", f)

check_sdf_matches(dg, "../../data/FormoseTestSetob.sdf")