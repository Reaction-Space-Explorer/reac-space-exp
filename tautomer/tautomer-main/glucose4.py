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
            >> rightPredicate[pred](repeat[3](inputRules))
	,ignoreRuleLabelTypes=True
        )


print("REACTION")

i= 0
reactants = {}
products = {}
rules=[]
for e in dg.edges: 
    reactants["reaction_" + str(i)] = [v.graph.smiles for v in e.sources]
    products["reaction_" + str(i)] =  [v.graph.smiles for v in e.targets]
    print([v.graph.smiles for v in e.sources],"reaction_"+str(i),[v.graph.smiles for v in e.targets])
    
    d = Derivations()
    d.rules = e.rules
    rules.append(d)
    print(d)
    i+=1

import pickle
dbfile = open('reactants', 'ab') 
pickle.dump(reactants, dbfile)
dbfile.close() 

dbfile = open('products', 'ab') 
pickle.dump(products, dbfile)
dbfile.close() 
  
  
  
print(reactants)




   







