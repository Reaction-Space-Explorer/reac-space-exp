import os
import time
#import compare_ms

include(os.path.abspath(os.path.join('..', 'rules/all.py')))
#include('clean_tautomers.py')

postChapter('Alkaline Glucose Degradation')

# starting molecule
#glucose = smiles('OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O', name='Glucose')
open_glucose = smiles("O=C[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)", "Open Chain Glucose")
water = smiles("O", name="Water")
# Forbidden substructures (don't add them to inputGraphs)
# Three and four membeered rings are unstable, C=C=C is forbidden
forbidden = [smiles('[C]1[C][C]1', name='cyclopropane', add=False), smiles('[C]1[C][C][C]1', name = 'cyclobutane', add=False),
            smiles('[C]1[C]O1', name='oxirane', add=False), smiles('[C]1[C][N]1',name='aziridine', add=False),
             smiles('[C]=[C]=[C]', name="Two double bonds", add=False)]


def pred(derivation):
    """
    Keyword arguments:
    d --- a derivation graph object
    """
    for g in derivation.right:
        # Allow masses only < 1200
        if g.exactMass >= 1200:
            return False
        for fb in forbidden:
            if fb.monomorphism(g) > 0:
                print(f"Found {fb} in {g}")
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

postSection('Final Network')
dg = DG(graphDatabase=inputGraphs)

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
        #subset, universe = clean_taut(dg, res)
        subset, universe = res.subset, res.universe
        #print('Subset size after removal:', len(subset))
        # This step replaces the previous subset (containing tautomers) with the cleaned subset
        #res = b.execute(addSubset(subset) >> addUniverse(universe))
        # now compare how many of these simulations were found in the MS data.
        #compare_ms.compare_sims([v.graph.smiles for v in dg.vertices], gen+1)
    print('Completed')

for e in dg.edges:
    for rule in e.rules:
        if "Benzilic" in rule.name:
            dg2 = DG(graphDatabase=inputGraphs)
            with dg2.build() as b:
                d = Derivations()
                sources = [source.graph for source in e.sources]
                targets = [target.graph for target in e.targets]
                d.left = sources
                d.rules = [rule]
                d.right = targets
                fake_edge = b.addDerivation(d)
            dg2.print()

print("Reaction rules that were used in this network generation were:")
rules_dict = {}
'''for e in dg.edges:
    for rule in e.rules:
        rules_dict[rule] = True
for rule in rules_dict.keys():
    print(rule)'''

# Print reactions of just the Benzilic Acid Rearrangement reaction alone
# dump smiles
'''with open("dump_smiles.txt", "w") as dump:
    for v in dg.vertices:
        dump.write(f'{v.graph.smiles}\n')'''

# Make a mass spectra (a histogram of the masses) of the molecules
#compare_ms.make_mass_spectra([v.graph.smiles for v in dg.vertices])
#dg.print()
'''postSection('Individual Vertices')
p = GraphPrinter()
p.simpleCarbons = True
p.withColour = True
p.collapseHydrogens = True
'''