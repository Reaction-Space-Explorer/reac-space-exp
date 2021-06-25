#This all6.py is a modification that seeks to generate the mass graphs to establish molecular mass evolution diagrams.
from matplotlib import pyplot
import math
include("common.py")



#Reaction to be inverted 
if True:
	include("sustract-proton.py")
	
##Check and Inversión	
canNotBeInvertedYet = []
loaded = list(inputRules)
for a in loaded:
	lt = a.labelType
	if lt is None:
		def checkLabel(l):
			if l == "*":
				print("Missing labelType in rule:", a.name)
				a.print()
				assert False
		for v in a.vertices:
			vl = v.left
			vr = v.right
			if not vl.isNull():
				checkLabel(vl.stringLabel)
			if not vr.isNull():
				checkLabel(vr.stringLabel)
		for e in a.edges:
			el = e.left
			er = e.right
			if not el.isNull():
				checkLabel(el.stringLabel)
			if not er.isNull():
				checkLabel(er.stringLabel)
	try:
		inv = a.makeInverse()
		if inv.isomorphism(a) == 0:
			inputRules.append(inv)
	except LogicError as e:
		print(a.name, "can not be inverted yet.")
		canNotBeInvertedYet.append(a)


#Reaction wihtout inversion
if True:

	include("OHRadO=M.py")
	include("doublebond-attack-Mrad.py")
	include("OHRadHM.py")
	include("termination.py")

##only Check 	
canNotBeInvertedYet = []
loaded = list(inputRules)
for a in loaded:
	lt = a.labelType
	if lt is None:
		def checkLabel(l):
			if l == "*":
				print("Missing labelType in rule:", a.name)
				a.print()
				assert False
		for v in a.vertices:
			vl = v.left
			vr = v.right
			if not vl.isNull():
				checkLabel(vl.stringLabel)
			if not vr.isNull():
				checkLabel(vr.stringLabel)
		for e in a.edges:
			el = e.left
			er = e.right
			if not el.isNull():
				checkLabel(el.stringLabel)
			if not er.isNull():
				checkLabel(er.stringLabel)
#	try:
#		inv = a.makeInverse()
#		if inv.isomorphism(a) == 0:
#			inputRules.append(inv)
#	except LogicError as e:
#		print(a.name, "can not be inverted yet.")
#		canNotBeInvertedYet.append(a)

#Define maximum number of heavy atoms and maximum number of radicals
def pred(derivation):
	for p in derivation.right:
		numAtoms = 0
		numRadicales = 0		
		for v in p.vertices:
			if v.stringLabel == "C":
				numAtoms += 1
			if v.stringLabel == "C.":
				numAtoms += 1			
			if v.stringLabel == "N":
				numAtoms += 1
			if v.stringLabel == "N.":
				numAtoms += 1			
			if v.stringLabel == "S":
				numAtoms += 1
			if v.stringLabel == "S.":
				numAtoms += 1	
			if v.stringLabel == "O":
				numAtoms += 1
			if v.stringLabel == "O.":
				numAtoms += 1					
				
		
			if v.stringLabel == "C.":
				numRadicales += 1
			if v.stringLabel == "N.":
				numRadicales += 1			
			if v.stringLabel == "S.":
				numRadicales += 1				
			if v.stringLabel == "O.":
				numRadicales += 1	
				
									
		if (numAtoms > 12 or numRadicales > 1) :
			return False
		
	return True
	
#Printing setup
p = DGPrinter()
p.withGraphName = False
#Read Etapa0.txt file for star the iteration
etapa=0
molecules = open('Etapa'+str(etapa)+'.txt','r')
lines = molecules.readlines()
print('lines'+str(lines))
namenumber=0
for line in lines:
	namenumber +=1
	line=line.rstrip('\n')
	mo = smiles(line, name="mol"+str(namenumber))
#	inputGraphs.append(mo)


for a in inputGraphs:
	print (a.smiles)
#Write Masas0.txt for the first Etapa0.txt
for g in inputGraphs:
	with open('Masas'+str(etapa)+'.txt', 'a') as the_file:
		the_file.write(str(g.exactMass)+'\n')



#Start iteration		
ciclo = True
dgdumslist=[]
while ciclo:
	print(etapa)
	#Plot molecular weights
	
	masas=[]
	molecules = open('Masas'+str(etapa)+'.txt','r')
	lines = molecules.readlines()
	for line in lines:
		masas.append(float(line))
	maxim=max(masas)
	minim=min(masas)
	binsfloat=maxim-minim
	binsinter=int(math.modf(binsfloat)[1]) +1 
	pyplot.hist(masas, bins=binsinter)
	pyplot.savefig('Densidad'+str(etapa))
	pyplot.clf()	
	
	#DG Graph generatino
	etapa+=1
	lenint= len(inputGraphs)
	print('inicialinputGraphs'+str(lenint))	
	#Ejecutar
	dg = DG(graphDatabase=inputGraphs)
	with dg.build() as b: 
		b.execute(
			addSubset(inputGraphs)
			>> 
			rightPredicate[
				pred
			]

			(
				repeat [1] (
			inputRules )
			),ignoreRuleLabelTypes=True
		)	
	print ('outputinputgraphs:'+str(len(inputGraphs)))
	#Print network
#	p.pushVertexVisible(lambda g, dg: g.smiles !='[H.]' and g.smiles !='[CH3.]' and g.smiles !='[HO.]' ) #Uncoment only when i need to erase somes species
	
	


	#EScribe la salida del prcesos ejecutado como smiles
	f= dg.dump()
	dg1=dgDump(inputGraphs,inputRules,f)
	for v in dg1.vertices:
#		v.graph.print()
#		print(v.graph.smiles)
		with open('Etapa'+str(etapa)+'.txt', 'a') as the_file:
			the_file.write(v.graph.smiles+'\n')
		with open('Masas'+str(etapa)+'.txt', 'a') as the_file2:
			the_file2.write(str(v.graph.exactMass)+'\n')
	
	
	
	""" get export files for importing into Neo4j by iterating across the hypergraph's nodes (molecules) and rels (edges) """
	dg_files_path = "Neo4j_Imports"
	generation_num = str(etapa)
	
	# write nodes txt file
	with open(dg_files_path + '/nodes/nodes_' + generation_num + '.txt', 'a') as nodes_file:
		for v in dg1.vertices:
			# node_id, smiles_str, exact_mass, node_label
			nodes_file.write(str(v.graph.id) + "," + v.graph.smiles + "," + str(v.graph.exactMass) + ",Molecule" "\n")
	
	# write relationships (rels) text file	
	with open(dg_files_path + "/rels/rels_" + generation_num + ".txt", 'a') as rels_file:
		for e in dg1.edges:
			print(e)
			for r in e.rules:
				print("\tRule:" + str(r))
			for s in e.sources:
				print("\tSource:" + str(s) + "," + str(s.graph.smiles))
			for t in e.targets:
				print("\tTarget:" + str(t) + "," + str(t.graph.smiles))
			
			# only one target per edge, but is collected in a list; same with reaction rule
			target = list(e.targets)[0]
			rule = list(e.rules)[0]
			
			# write line for each target molecule
			for source in e.sources:
				source_smiles = str(source.graph.smiles)
				target_smiles = str(target.graph.smiles)
				reaction_rule = str(rule)
				# edge_id, source_smiles_str, target_smiles_str, reaction_rule_str
				rels_file.write(str(e.id) + "," + source_smiles + "," + target_smiles + "," + reaction_rule + "\n")
	
	""" end Neo4j export block """
	
	
	

	inputGraphs=[] # We need this to avoid isomorphism problems between new outputs
	
	#Leer nuevamente las moleculas en la etapa
	namenumber=0
	molecules = open('Etapa'+str(etapa)+'.txt','r')
	lines = molecules.readlines()
	for line in lines:
		namenumber +=1
		line=line.strip()
		mo = smiles(line, name="mol"+str(namenumber))
		inputGraphs.append(mo)
	
	dg.print(p)
	
	#Append name of the dgdump files	
	dgdumslist.append(f)

	#copy dgdump files to my folder
	from shutil import copyfile
	copyfile(f, "dgdumfile-"+str(etapa))

	#Condition to finish loop	
	if (len(inputGraphs) ==  lenint):
		ciclo=False

print("dump list of file is:",dgdumslist)
for a in inputRules:
	a.print()

