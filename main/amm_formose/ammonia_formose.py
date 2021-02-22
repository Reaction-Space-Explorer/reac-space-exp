# NH3 + HCHO + HOCH2CHO

# use formaldehyde as a component in Cannizarro, as opposed to the rule for glucose Cannizarro
with_formaldehyde=True

include("../main.py")

postChapter("Ammonia Formose Reaction")

ammonia = smiles("N", name="Ammonia")
formaldehyde = smiles("C=O", name="Formaldehyde")
glycoladehyde = smiles("OCC=O", name="Glycolaldehyde")
water = smiles("O", name="Water")

generations = 4

dg = DG(graphDatabase=inputGraphs,
	labelSettings=LabelSettings(LabelType.Term, LabelRelation.Specialisation))

subset = inputGraphs
universe = []
# dump initial reactants as part of "G0"
write_gen_output(subset, generation=0, reaction_name="formose_amm")

with dg.build() as b:
	for gen in range(generations):
		start_time = time.time()
		print(f"Starting round {gen+1}")
		res = b.execute(addSubset(subset) >> addUniverse(universe) >> strat, verbosity=2)
		end_time = time.time()
		print(f"Took {end_time - start_time} seconds to complete round {gen+1}")
		print(f'Products in generation {gen+1}:', len(res.subset))
		subset, universe = res.subset, res.universe
		export_to_neo4j(dg_obj = dg, generation_num = gen+1)
		write_gen_output(subset, gen+1, reaction_name="formose_amm")
	print('Completed')

# Dump the dg so it can be loaded again quickly without having to generate it from scratch.
f = dg.dump()
print("Dump file: ", f)
