import os
import time
import compare_ms
from rdkit.Chem import SDMolSupplier, MolToSmiles
from mod_to_neo4j_exporter import export_to_neo4j

include(os.path.abspath(os.path.join('..', 'rules/all.py')))
#include('clean_tautomers.py')

postChapter('Alkaline Glucose Degradation')

# starting molecule
#glucose = smiles('OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O', name='Glucose')
open_glucose = smiles("O=C[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)", "Open Chain Glucose")
water = smiles("O", name="Water")
# Forbidden substructures (don't add them to inputGraphs)
# Three and four membeered rings are unstable, any atom with two double bonds is forbidden
# Note: the latter means that even O=C=O (carbon dioxide) will also be forbidden
forbidden = [smiles('[*]1[*][*]1', name="Three membered ring", add=False),
             smiles("[*]1[*][*][*]1", name="Four membered ring", add=False),
             smiles('[*]=[*]=[*]', name="Two double bonds", add=False)]

def pred(derivation):
    """
    Keyword arguments:
    d --- a derivation graph object
    """
    for g in derivation.right:
        # Allow masses only < 1200
        if g.exactMass >= 500:
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
generations = 2

#postSection('Final Network')
dg = DG(graphDatabase=inputGraphs,
    labelSettings=LabelSettings(LabelType.Term, LabelRelation.Specialisation))

'''dg = dgDump(inputGraphs, inputRules, "000_DG.dg")
print("Finished loading from dump file")'''

subset = inputGraphs
universe = []
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
        # now compare how many of these simulations were found in the MS data.
        compare_ms.compare_sims([v.graph.smiles for v in dg.vertices], gen+1)
        export_to_neo4j(dg_obj = dg, generation_num = gen)
    print('Completed')

# compare structures with what the Y&M paper has
sdfile = SDMolSupplier(os.path.join("..", "data/NewAlkalineHydrolysisStructures.sdf"))

matching_structs = []
not_matching = []

print("Checking for matches with Y&M's structures")

for mol in sdfile:
    mol_graph = smiles(MolToSmiles(mol), add=False)
    for g in universe:
        if g.isomorphism(mol_graph) > 0:
            matching_structs.append(mol_graph)
            print("Structure {0} of the SDF found in the network!".format(mol_graph.smiles))
        else:
            not_matching.append(mol_graph)

print(f"{100* len(matching_structs)/len(sdfile)}% structures in the SDF are in the reaction network.")

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
compare_ms.make_mass_spectra([v.graph.smiles for v in dg.vertices])

# Compare structures with 

# Print reactions of just the one reaction alone, in separate DGs
to_print = ['Cannizarro', 'Knoevenagel, C, C(=O)A, OR, H |, HC', 'Ester Formation Hydrolysis Exchange',
            'Knoevenagel, C, C(=O)A, OR, H |, CC', 'Ring Closure']

dgprint = DGPrinter()
dgprint.withRuleName = True
dgprint.withShortcutEdges = True

for item_to_print in to_print:
    count = 0
    postSection(f"{item_to_print} reactions")
    for e in dg.edges:
        # Don't print more than 35 of any category
        if count > 35:
            pass
        #    break
        else:
            for rule in e.rules:
                if item_to_print in rule.name:
                    count += 1
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
                    dg2.print(dgprint)
# dump smiles
'''with open("dump_smiles.txt", "w") as dump:
    for v in dg.vertices:
        dump.write(f'{v.graph.smiles}\n')'''

#dg.print()
'''postSection('Individual Vertices')
p = GraphPrinter()
p.simpleCarbons = True
p.withColour = True
p.collapseHydrogens = True
'''