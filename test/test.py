g1 = smiles('C=C=C=C=C', name="Test", add=False)
#g1.print()
'''dumb_rule = ruleGMLString("""rule [
    ruleID "Dumb Rule"
    left [
        edge [ source 0 target 1 label "="]
        ]
    context [
        node [ id 0 label "C"]
        node [ id 1 label "C"]
    ]
    right [
        edge [ source 0 target 1 label "-"]
    ]
]""")
dumb_rule.print()'''

#include(os.path.join("..", "rules/benzilicAcidRearrangement.py"))
#include(os.path.join("..", "rules/aldolCondensationOneStep.py"))
#include(os.path.join("..", "rules/aldolCondensation.py"))
include(os.path.join("..", "rules/elimination1.py"))
include(os.path.join("..", "rules/elimination2.py"))
include(os.path.join("..", "rules/retroAldol.py"))
include(os.path.join("..", "rules/ketoEnolisation.py"))
include(os.path.join("..", "rules/hydration.py"))
ketoEnolisation[0].makeInverse()
water = smiles("O", name="Water")
#test_molecule = smiles("O=C(C(=O)c1ccccc1)c1ccccc1", name = "Test Aromatic Compound")
#cyclohexandione = smiles("O=C1CCCCC1=O", name="Cyclohexan-1,2-dione")
#straight_chain = smiles("C(=O)C(=O)", name="Pentan-2,3-dione")
#acetaldehyde = smiles("CC=O", name="Acetaldehyde")
#ethanamide = smiles("NC(=O)C", name = "Ethanamide")
#mercaptan = smiles("SC(=O)C", name = "Mercaptoethanone")
#formaldehyde = smiles("C=O", name="Formaldehyde")
aldol = smiles("O=CCC(O)C", name='Aldol')
p = GraphPrinter()
p.collapseHydrogens = True
p.simpleCarbons = True
p.withColour = True

#test_molecule.print(p)

forbidden = graphGMLString("""graph[
    node [ id 0 label "*" ]
    node [ id 1 label "*" ]
    node [ id 2 label "*"]
    edge [ source 0 target 1 label "=" ]
    edge [ source 1 target 2 label "=" ]
]
""", name="Two double bonds", add=False)

# Test if double bonds can 
def pred(derivation):
    for g in derivation.right:
        if forbidden.monomorphism(g, 
            labelSettings=LabelSettings(LabelType.Term, LabelRelation.Specialisation)) > 0:
            print("Found double bonds!")
            return False
    return True

for rule in inputRules:
    rule.print()

generations = 3

strat = (
    addSubset(inputGraphs)
    >> rightPredicate[pred](repeat[generations](inputRules))
)


dg = DG(graphDatabase=inputGraphs, labelSettings=LabelSettings(LabelType.Term,LabelRelation.Specialisation))
with dg.build() as b:
    b.execute(strat, verbosity=8)

dg.print()

postSection("Elimination")
count = 0
for e in dg.edges:
    if count > 100:
        break
    else:
        for rule in e.rules:
            if "Elimination +" in rule.name:
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
                dg2.print()

postSection("Individual Molecules in the Network")
for g in dg.vertices:
    g.graph.print(p)

'''
def dumpDerivation(graph, filename):
    # graph is a DG object
    with open(filename, 'w') as file:
        file.write("# (Format) Vertex: <id>, <graphDFS>, <inEdges>, <outEdges> \n")
        file.write("# Edge: <id>, <ruleName> ")
        for v in graph.vertices:
            file.write(f'Vertex:{v.id}\t{v.graph.graphDFS}\t{[e.id for e in v.inEdges]}\t{[e.id for e in v.outEdges]}')
            file.write('\n')
        for e in graph.edges:
            # save edge id and rule names
            file.write(f'Edge:\t{e.id}\t{[r.name for r in e.rules]}')
            file.write('\n')
    print('Dump file: ', filename)

def loadDump(filename):
    try:
        with open(filename) as file:
            for line in file.readlines():
                if line.startswith('#'):
                    # Ignore comments
                    pass
                lineElements = line.split(sep='\t')
                if lineElements[0] == 'Edge:':
                    pass
    except:
        print(f'Error writing {filename}')
'''
'''some_rule = ruleGMLString("""rule[
    left[
        node [ id 0 label "C"]
    ]
]
""")
some_rule.print()
dg = DG(graphDatabase=inputGraphs)
with dg.build() as b:
    res = b.execute(strat)
    subset = res.subset
    universe = res.universe
    print(universe)
    d = Derivations()
    d.left = [universe[0]]
    d.rules = [some_rule]
    d.right = [universe[1]]
    new_edge = b.addDerivation(d)
    #new_edge.print()
    for i in range(generations-1):
        subset = res.subset
        for item in subset:
            print(type(item))
        universe = res.universe
        b.execute(addSubset(subset) >> addUniverse(universe) >> strat)
dg.print()'''

#dumpDerivation(dg, 'test_dump.txt')