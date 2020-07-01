import os
from rdkit.Chem import MolFromSmiles, MolToSmiles

include(os.path.abspath(os.path.join('..', 'rules/all.py')))
include('clean_tautomers.py')

postChapter('Alkaline Glucose Degradation')

# starting molecule
glucose = smiles('OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O', name='Glucose')

# Forbidden substructures
# Three and four membeered rings are unstable, C=C=C is forbidden
forbidden = [smiles('[C]1[C][C]1', name='cyclopropane'), smiles('[C]1[C][C][C]1', name = 'cyclobutane'),
            smiles('[C]1[C]O1', name='oxirane'), smiles('[C]1[C][N]1',name='aziridine'), smiles('[C]=[C]=[C]', name="Two double bonds")]

# make sure these don't get passed as an input
for fb in forbidden:
    inputGraphs.remove(fb)

def pred(derivation):
    """
    Keyword arguments:
    d --- a derivation graph object
    """
    for g in derivation.right:
        if g.exactMass >= 500:
            return False
        for fb in forbidden:
            if fb.monomorphism(g) > 0:
                print(f"Found {fb} in {g}")
                return False
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
with dg.build() as b:
    res = b.execute(strat, verbosity=2, ignoreRuleLabelTypes=True)
    subset, universe = clean_taut(dg, res)
    for gen in range(generations-1):
        res = b.execute(addSubset(subset) >> addUniverse(universe) >> strat,
                            verbosity=2, ignoreRuleLabelTypes=True)
        print('Original subset size:', len(res.subset))
        subset, universe = clean_taut(dg, res)
        print('Subset size after removal:', len(subset))
        # This step replaces the previous subset (containing tautomers) with the cleaned subset
        res = b.execute(addSubset(subset) >> addUniverse(universe))
    print('Completed')
#dg.print()
postSection('Individual Vertices')
p = GraphPrinter()
p.simpleCarbons = True
p.withColour = True
p.collapseHydrogens = True
'''
for v in dg.vertices:
    v.graph.print(p)
postSection('Individual Edges')
for e in dg.edges:
    e.print()'''