#This all7.py is a modification that seeks to avoid forbidden chains and also include Jessica's Neo4J code.
from matplotlib import pyplot
import math
################### Commun functions ################################################
#include("common.py")

################### Clear previous output ###########################################
import os
import sys
import glob
import os.path

##remove Mass and list of molecules
list_of_files = glob.glob('data/*.txt') 
for file_name in list_of_files:
	try:
		f = open(file_name)
		os.remove(file_name)
	except IOError:
	    print("File not accessible")

##remove plots
list_of_files = glob.glob('data/*.png') 
for file_name in list_of_files:
	try:
		f = open(file_name)
		os.remove(file_name)
	except IOError:
	    print("File not accessible")

##remove Neo4j output
list_of_files = glob.glob('Neo4j_Imports/nodes/*.txt') 
for file_name in list_of_files:
	try:
		f = open(file_name)
		os.remove(file_name)
	except IOError:
	    print("File not accessible")

list_of_files = glob.glob('Neo4j_Imports/rels/*.txt') 
for file_name in list_of_files:
	try:
		f = open(file_name)
		os.remove(file_name)
	except IOError:
	    print("File not accessible")


################### forbidden molecules #############################################
include("forbidden.py")

################ Reaction rules #####################################################
if True:
	include("rules/ozone-rad1.py")
	include("rules/ozone-rad2.py")
	include("rules/ozone-rad3.py")
	include("rules/h-abstr-rad1.py")
	include("rules/h-abstr-rad2.py")
	include("rules/h-abstr-rad3.py")
	include("rules/oxygen-rad1.py")
	include("rules/oxygen-rad2.py")
	include("rules/termo-oxygen-rad1.py")
	include("rules/termo-oxygen-rad2.py")
	include("rules/termo-oxygen-rad3.py")
	include("rules/termo-hy-rad1.py")
	include("rules/termo-hy-rad2.py")
	include("rules/foto-1.py")	
	include("rules/foto-2.py")
	include("rules/foto-3.py")
	include("rules/termination-rad1.py")

#################  inputGraphs redefinition ##########################################

inputGraphs = [] # In order to avoid the automatic input of molecules in the glabal list inputGraphs

################### Termolecular reaction and radiation #############################

if True: 
	termolecular = smiles("[Md]", name="termolecular")
	inputGraphs.append(termolecular)
if True:	
	radiation = smiles("[Hf]", name = "radiation")
	inputGraphs.append(radiation)

#################  Condition for dg generation #######################################
#Define maximum number of heavy atoms, maximum number of radicals and forbidden structures
def pred(derivation):

	for p in derivation.right:
		numAtoms = 0
		numRadicales = 0
		numOx = 0
		 

		numchainmatch = 0		
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



			
		for forb in forbidden:
			mor1 = forb.monomorphism(p, maxNumMatches=2**30,labelSettings=LabelSettings(LabelType.Term, LabelRelation.Specialisation))
			isom1 = forb.isomorphism(forb, maxNumMatches=2**30,labelSettings=LabelSettings(LabelType.Term, LabelRelation.Specialisation))
			numchainmatch += mor1/isom1

		
		
									
		if (numAtoms > 10 or numRadicales > 1 or numchainmatch > 0) :
			return False	
	return True
	
################### Printing network setup #########################################
p = DGPrinter()
p.withGraphName = False
p.withRuleId = True

#p.pushVertexVisible(lambda g, dg: g.smiles !="[Md]") 

##################################################################################
################### Beford start network generation ###############################
###################################################################################

###################Read Etapa0.txt file ###########################################
#for star the iteration and upload the first molecules to inputGraphs list
etapa=0
molecules = open('Etapa'+str(etapa)+'.txt','r')
lines = molecules.readlines()
print('lines'+str(lines))
for line in lines:
	line=line.rstrip('\n')
	mo = smiles(line)
	inputGraphs.append(mo)
for a in inputGraphs:
	print (a.smiles)


#################### Write Masas0.txt ################################################
#for the first Etapa0.txt
for g in inputGraphs:
	with open('data/Masas'+str(etapa)+'.txt', 'a') as the_file:
		if (g.smiles != "[Hf]" and g.smiles != "[Md]"):
			the_file.write(str(g.exactMass)+'\n') 


#######################################################################################
#################### "" Start iteration ""	#########################################
########################################################################################	
ciclo = True
dgdumslist=[]
postSection("DG graphs")
while ciclo:
	print('Step: ', etapa)
	############### Plot molecular weights ##########################################
	
	masas=[]
	molecules = open('data/Masas'+str(etapa)+'.txt','r')
	lines = molecules.readlines()
	for line in lines:
		masas.append(float(line))
	maxim=max(masas)
	minim=min(masas)
	binsfloat=maxim-minim
	binsinter=int(math.modf(binsfloat)[1]) +1 
	pyplot.hist(masas, bins=binsinter)
	pyplot.savefig('data/Densidad'+str(etapa))
	pyplot.clf()	
	
	############### DG Graph generation  ################################################
	etapa+=1
	lenint= len(inputGraphs)
	#Ejecutar
	dg = DG(graphDatabase=inputGraphs,labelSettings=LabelSettings(LabelType.Term,LabelRelation.Specialisation))
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
			)
		)

	for v in inputGraphs:
		print (v.smiles)
	
	dg.print(p)
	################## Write the output of the step #######################################
	
	f = dg.dump()
	dg1=dgDump(inputGraphs,inputRules,f)
	for v in dg1.vertices:
#		v.graph.print()
#		print(v.graph.smiles)
		with open('data/Etapa'+str(etapa)+'.txt', 'a') as the_file:
			the_file.write(v.graph.smiles+'\n')
		with open('data/Masas'+str(etapa)+'.txt', 'a') as the_file2:
			if (v.graph.smiles != "[Hf]" and v.graph.smiles != "[Md]"):
				the_file2.write(str(v.graph.exactMass)+'\n')

	###################### Neo4J #############################################################
	""" get export files for importing into Neo4j by iterating across the hypergraph's nodes (molecules) and rels (edges) """
	dg_files_path = "Neo4j_Imports"
	generation_num = str(etapa)
	
	# write nodes txt file
	with open(dg_files_path + '/nodes/nodes_' + generation_num + '.txt', 'a') as nodes_file:
		for v in dg1.vertices:
			#print(v)
			# node_id, smiles_str, exact_mass, node_label #Romulo: termolecular species Md and radiation Hf has ficticious mass
			if (v.graph.smiles != "[Hf]" and v.graph.smiles != "[Md]"):
				nodes_file.write(str(v.id) + "," + v.graph.smiles + "," + str(v.graph.exactMass) + ",Molecule" "\n")
			
			if (v.graph.smiles == "[Hf]" ): 
				nodes_file.write(str(v.id) + "," + v.graph.smiles + "," + "Without mass" + ",Radiation" "\n")
			
			if (v.graph.smiles == "[Md]" ): 
				nodes_file.write(str(v.id) + "," + v.graph.smiles + "," + "Mass of the most abundat" + ",Termolecular" "\n")
					
	# write relationships (rels) text file	
	with open(dg_files_path + "/rels/rels_" + generation_num + ".txt", 'a') as rels_file:
		print("#################### Inicio de dg.edges#####################")
		print (str(dg.edges))
		print("#################### Fin de dg.edges#####################")
		for e in dg1.edges:
			print("edge id: " + str(e.id))
			print("For neo4j:")
			print(e)
			for r in e.rules: # e.rules has only one element
				print("\tRule:" + str(r)+ "," +str(r.id))
			for s in e.sources:
				print("\tSource:" + str(s) + "," + str(s.graph.smiles)+"," + str(s.id))
			for t in e.targets:
				print("\tTarget:" + str(t) + "," + str(t.graph.smiles)+"," + str(t.id))
			
			
			# only one target per edge, but is collected in a list; same with reaction rule
			#target = list(e.targets)[0]
			rule = list(e.rules)[0]
			#print(str(e.rules))
			
			for target in e.targets:
			# write line for each target molecule
				for source in e.sources:
					source_smiles = str(source.graph.smiles)
					target_smiles = str(target.graph.smiles)
					reaction_rule = str(rule)
					# edge_id, source_smiles_str, target_smiles_str, reaction_rule_str
					rels_file.write(str(e.id) + "," + source_smiles + "," + target_smiles + "," + reaction_rule + "\n")

	""" end Neo4j export block """
	
	
	
	##########################################################################################
	
	
	
	inputGraphs=[] # We need this to avoid isomorphism problems between new outputs
	
	
	
	
	#################  Read the output of the step  ######################################## 
	molecules = open('data/Etapa'+str(etapa)+'.txt','r')
	lines = molecules.readlines()
	for line in lines:
		line=line.strip()
		mo = smiles(line)
		inputGraphs.append(mo)
	
	
	

	################# Condition to finish loop  ################################################# 	
	if (len(inputGraphs) ==  lenint):
		ciclo=False

#######################################################################################
#################### "" End iteration ""	#########################################
########################################################################################


print("dump list of file is:",dgdumslist)
postSection("Rules")
for a in inputRules:
	a.print()

