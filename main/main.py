import os

include(os.path.abspath(os.path.join('..', 'rules/all.py')))
include('clean_tautomers.py')

postChapter('Alkaline Glucose Degradation')

# starting molecule
glucose = smiles('OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O', name='Glucose')

# Forbidden substructures
# Three and four membeered rings are unstable
oxirane = graphGMLString("""graph [
   node [ id 0 label "C" ]
   node [ id 1 label "C" ]
   node [ id 2 label "O" ]
   edge [ source 0 target 1 label "-"]
   edge [ source 0 target 2 label "-"]
   edge [ source 1 target 2 label "-"]
]""", name='Oxirane')

aziridine = graphGMLString("""graph [
   node [ id 0 label "C" ]
   node [ id 1 label "C" ]
   node [ id 2 label "N" ]
   edge [ source 0 target 1 label "-"]
   edge [ source 0 target 2 label "-"]
   edge [ source 1 target 2 label "-"]
]""", name='aziridine')
forbidden = [smiles('C1CC1', name='cyclopropane'), smiles('C1CCC1', name = 'cyclobutane'),
            oxirane, aziridine]

# make sure these don't get passed as an input
for fb in forbidden:
    inputGraphs.remove(fb)
strat = (
    addSubset(inputGraphs)
    >> rightPredicate[
        lambda derivation: all (g.exactMass <= 500
         and (g.monomorphism(fb) == 0 for fb in forbidden) for g in derivation.right)
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
'''for fb in forbidden:
    fb.print()
for v in dg.vertices:
    for item in forbidden:
        if v.graph.monomorphism(item) == 1:
            print(f'Found substructure {item.name} in {v.graph.name}')'''
'''for v in dg.vertices:
    v.graph.print(p)
postSection('Individual Edges')
for e in dg.edges:
    e.print()'''