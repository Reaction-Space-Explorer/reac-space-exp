# There's no formaldehyde at the beginning, turn the Cannizarro rule for formaldehyde off
with_formaldehyde = False

include("../main.py")
# including these files using MOD's include() so that MOD's functions are callable in them
include("../mod_to_neo4j_exporter.py")
#include('clean_tautomers.py')

postChapter('Alkaline Glucose Degradation')

open_glucose = smiles("O=CC(O)C(O)C(O)C(O)C(O)", "Open Chain Glucose")
water = smiles("O", name="Water")

# Number of generations we want to perform
generations = 3

print("Loading .DG file")
start_t = time.time()

dg = DG.load(inputGraphs, inputRules, "glucose_5g_21jan.dg")
print(f"Finished loading from dump file, took {time.time() - start_t} seconds")

print(f'Number of edges in DG: {len(list(dg.edges))}')
rule_count = [len(list(e.rules)) for e in dg.edges]
print(f'First 20 elements of rule_count list: {rule_count[:10]}')
print(f'Total: {sum(rule_count)}')

count_rules_by_gen(dg, 'glucose_degradation_output_10mar.txt')

'''dg_2 = DG.load(dg.graphDatabaseexit, inputRules, "glucose_5g_10mar.dg")
print(f"Finished loading from dump file, took {time.time() - start_t} seconds")
'''

#diffDGs(dg, dg_2)

# TODO: make this a more general method, which can be called in other reactions.
'''print('Exporting reactions list for thermo calculations.')
with open('reactions_list_smi.txt', 'w') as reac_list:
	for edge in dg.edges:
		#print(f"Edge: {edge.id}")
		source_smiles = [s.graph.smiles for s in edge.sources]
		#print('Sources:', source_smiles)
		target_smiles = [t.graph.smiles for t in edge.targets]
		#print('Targets', target_smiles)
		reaction_str = f"{' + '.join(source_smiles)} -> {' + '.join(target_smiles)}"
		reac_list.write(f'{reaction_str}\n')'''


'''dg = DG(graphDatabase=inputGraphs,
	labelSettings=LabelSettings(LabelType.Term, LabelRelation.Specialisation))

subset = inputGraphs
universe = []

# dump initial reactants as part of "G0"
write_gen_output(subset, generation=0, reaction_name="glucose_degradation")

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
		#export_to_neo4j(dg_obj = dg, generation_num = gen+1)
		write_gen_output(subset, gen+1, reaction_name="glucose_degradation")
	print('Completed')'''

# Dump the dg so it can be loaded again quickly without having to generate it from scratch.
f = dg.dump()
print("Dump file: ", f)


'''check_sdf_matches(dg, "../../data/NewAlkalineHydrolysisStructures.sdf")

enol = smiles("[C]=[C]O", name="enol substruct", add=False)
diol = smiles("O[C]O", name="diol substruct", add=False)
find_substruct_producer(dg, enol, print_rule=True)
find_substruct_producer(dg, diol, print_rule=True)'''