# time: 1:10:00 at video
include("../hcn/rules/all.py") # mod way of importing files

#postSection("Input Rules") #CREATING A SECTION FOR RULES
#for a in inputRules:       # PRINTING RULES
#    a.print()



# INPUT MOLECULES (GRAPHS)
postSection("Input MOLECULES") 

glucose = smiles('OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O', name='Glucose')
water  = smiles('O', name ='water') 
# PRINT ALL inputGraphs (INPUT MOLECULES)
postSection("Input Molecules") 
for a in inputGraphs:
    a.print()

# add all input graphs to Derivation graph
dg = DG(graphDatabase=inputGraphs)

strat = (
    #addSubset()
    addSubset(glucose,water)
    >> rightPredicate[lambda derivation: all (g.exactMass <= 200 for g in derivation.right)](repeat[2](inputRules))
)

dg.build().execute(strat, ignoreRuleLabelTypes=True)
#dg.dump()
dg.print()
