# time: 1:10:00 at video
include("../hcn/rules/all.py") # mod way of importing files
# INPUT MOLECULES (GRAPHS)
postSection("Input MOLECULES") 

glucose = smiles('OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O', name='Glucose')
water  = smiles('O', name ='water') 
# PRINT ALL inputGraphs (INPUT MOLECULES)
postSection("Input Molecules") 
initMolecules = 0
for a in inputGraphs:
    initMolecules += 1
    a.print()

# add all input graphs to Derivation graph
dg = DG(graphDatabase=inputGraphs)

# right predicate
def pred(derivation):
    for p in derivation.right:
        numCs = 0
        for v in p.vertices:
            if v.stringLabel == "C":
                numCs += 1
        if numCs > 6:
            return False

        return True 
# APPLY DERIVATIONS
with dg.build() as b:
    b.execute(
            addSubset(glucose, water)
            >> rightPredicate[pred](repeat[1](inputRules))
	,ignoreRuleLabelTypes=True
        )

# SAVE NETWORK
#f = dg.dump()
#print("DUBED FILE SAVED AS:", f)

#pints network
dg.print()

# A pesonalized summary for me
print("\nMY OWN SUMMARY:")
totalMolecules = dg.numVertices
outMolecules = abs(initMolecules - totalMolecules)

print("TOTAL MOLECULES:",totalMolecules)
print("INITIAL MOLECULES:",initMolecules)
print("MOLECULES GENERATED:",outMolecules,"\n")

#print molecules in a txt file 
count = 0
with open("out_mol.txt",'w') as w:
    for v in dg.vertices:
        count += 1
        w.write('"')
        w.write(v.graph.smiles)
        w.write('"')
        if count < totalMolecules:
            w.write(",")
        print(v.graph.smiles)
