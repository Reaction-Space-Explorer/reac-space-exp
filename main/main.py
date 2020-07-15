import os
import time
import compare_ms
from rdkit.Chem import SDMolSupplier, MolToSmiles

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
forbidden = [smiles('[C]1[C][C]1', name='cyclopropane', add=False), smiles('[C]1[C][C][C]1', name = 'cyclobutane', add=False),
            smiles('[C]1[C]O1', name='oxirane', add=False), smiles('[C]1[C][N]1',name='aziridine', add=False),
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
                #print(f"Found {fb} in {g}")
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
dg = DG(graphDatabase=inputGraphs)

'''dg = dgDump(inputGraphs, inputRules, "000_DG.dg")
print("Finished loading from dump file")'''

subset = inputGraphs
universe = []
with dg.build() as b:
    for gen in range(generations):
        start_time = time.time()
        print(f"Starting round {gen+1}")
        res = b.execute(addSubset(subset) >> addUniverse(universe) >> strat,
                            verbosity=2, ignoreRuleLabelTypes=True)
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
print("Aldol condensation reaction count: ", rules_count.count("Aldol Condensation"))

print("Rules used: {0}".format(dict({rule:True for rule in rules_count}).keys()))
# Make a mass spectra (a histogram of the masses) of the molecules
compare_ms.make_mass_spectra([v.graph.smiles for v in dg.vertices])

# Compare structures with 

# Print reactions of just the one reaction alone, in separate DGs
'''count = 0
for e in dg.edges:
    if count > 100:
        pass
    #    break
    else:
        for rule in e.rules:
            if "Aldol Condensation" == rule.name:
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
                dg2.print()'''
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