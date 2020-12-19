with_formaldehyde = True
 
include("../main.py")
# including these files using MOD's include() so that MOD's functions are callable in them
include("../mod_to_neo4j_exporter.py")
#include('clean_tautomers.py')
 
postChapter('Pyruvic Acid')
 
pyruvic = smiles("CC(=O)C(=O)O", "Pyruvic Acid")
water = smiles("O", name="Water")# Number of generations we want to perform
generations = 5
 
dg = DG(graphDatabase=inputGraphs,
    labelSettings=LabelSettings(LabelType.Term, LabelRelation.Specialisation))
 
subset = inputGraphs
universe = []

# write gen 0 (starting material in output .txt file)
write_gen_output(subset, generation=0, reaction_name="pyruvic")

with dg.build() as b:
    for gen in range(generations):
        start_time = time.time()
        print(f"Starting round {gen+1}")
        res = b.execute(addSubset(subset) >> addUniverse(universe) >> strat, verbosity=8)
        end_time = time.time()
        print(f"Took {end_time - start_time} seconds to complete round {gen+1}")
        print(f'Products in gen {gen+1}:', len(res.subset))
        subset, universe = res.subset, res.universe
        export_to_neo4j(dg_obj = dg, generation_num = gen)
        write_gen_output(subset, gen+1, reaction_name="pyruvic")
    print('Completed')
 
# Dump the dg so it can be loaded again quickly without having to generate it from scratch.
f = dg.dump()
print("Dump file: ", f)

with open('pyruvic_4gens.txt', 'w') as wr:
    for v in dg.vertices:
        wr.write(f'{v.graph.smiles}\n')