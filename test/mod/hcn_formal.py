include("rules/all.py")

print('Hello World!')

# Make some headers in the summary:
postChapter("HCN-HCHO")

postSection("Starting Molecules")

hcn = smiles('C#N', name='HCN')
formaldehyde = smiles('C=O', name='Formaldehyde')
hcn.print()
formaldehyde.print()
glucose = smiles('OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O', name='Glucose')
glucose.print()
glycolaldehyde = smiles('OCC=O', name = 'Glycolaldehyde')
postSection("Reactions")

strat = (
    addSubset(formaldehyde)
    >> addSubset(glucose)
    >> rightPredicate[
        lambda derivation: all (g.exactMass <= 500 for g in derivation.right)
    ]
    (
        repeat[3](inputRules))
)
dg = DG(graphDatabase=inputGraphs)
dg.build().execute(strat, ignoreRuleLabelTypes=True)
#dg.dump()
#dg.print()


postSection("Resulting Molecules")
p = GraphPrinter()
p.edgesAsBonds = True
p.withColour = True
p.collapseHydrogens = True
p.simpleCarbons = True
output_file = open('generated_mols.txt','w')
output_file.write('# Note: this contains the input molecules as well \n')
for v in dg.vertices:
    output_file.write(f'{v.graph.smiles}\n')