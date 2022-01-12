from matplotlib import pyplot
import pandas as pd
import sys
import json
import os
import os.path
import glob
import numpy as np 
import struct
import csv
import ast
from random import sample
import math
from pathlib import Path
import shutil

# This program transform the rel_x.txt data and get a format necesary for gephi and also return a DG energy of reaction dictionary in case we have equilibrator data, else return a empty dictionary. Also return a dictionary for reactants and products in each raction node.

def get_data_pd(data_dir,file_name,spontan): 
	
	out_rel_name = data_dir+file_name[:-4]+"_gephi.csv"
	out_types_name = data_dir+file_name[:-4]+"_types_gephi.csv"
	file_name = data_dir+file_name

	### Erase previous outputs
	filed = Path(out_rel_name)
	#print(filed)
	#print(filed.is_file())
	if filed.is_file():
		question = input("Do you want to erase the previos files "+out_rel_name+ "," + out_types_name + "(yes/not):")
		if question == "yes" or question == "Yes" or question == "YES":
			os.remove(out_rel_name)
			os.remove(out_types_name)
		else:
			print("Finishing the script")
			exit()	
	### 
	molecules = open(file_name,'r')
	lines = molecules.readlines()
	### building rel_gephi output
	with open(out_rel_name, 'a') as reaction_file:
		reaction_file.write( "source"+"\t"+"target"+"\t"+"reac_name"+"\t"+"counter"+"\t"+"DG"+"\n")
		i =0
		for line in lines:
			left = []
			right = []

			line=line.rstrip('\n')
			a=line.split('\t')

			if a[2] =="1":
				left.append(a[0])
				right.append(a[1])
			if a[2] =="-1":
				left.append(a[1])
				right.append(a[0])
			i += 1
			if left[0][:-2].isnumeric():
				#if not left[0][:-2] in spontan:
				if not left[0] in spontan:
					#print("Error: the reaction: "+ left[0][:-2] + " Isn't in spontan, set correctly this ")
					print("Error: the reaction: "+ left[0] + " Isn't in spontan, set correctly this ")
					exit()
				#fre_energy = spontan[left[0][:-2]]				
				fre_energy = spontan[left[0]]
				
			if right[0][:-2].isnumeric():
				#if not right[0][:-2] in spontan:
				if not right[0] in spontan:
					#print("Error: the reaction: "+ right[0][:-2] + " Isn't in spontan, set correctly this ")
					print("Error: the reaction: "+ right[0] + " Isn't in spontan, set correctly this ")
					exit()	
				#fre_energy = spontan[right[0][:-2]]		
				fre_energy = spontan[right[0]]
			reaction_file.write( left[0]+"\t"+right[0]+"\t"+a[3]+"\t"+str(i)+"\t"+str(fre_energy)+"\n")

			
	### Building types gephi output
	reac_center = []
	molecu = []
	for line in lines:
		line=line.rstrip('\n')
		a=line.split('\t')
		reac_center.append(a[0])
		molecu.append(a[1])
	reac_center=list(set(reac_center))
	molecu=list(set(molecu))

	with open(out_types_name,'a') as reaction_file:
		reaction_file.write( "Id"+"\t"+"type"+"\n")
		for mol in molecu:
			reaction_file.write(mol+"\t"+"1"+"\n")	
		for reac in reac_center:
			reaction_file.write(reac+"\t"+"2"+"\n")
	print("Files created:")
	print(out_rel_name)
	print(out_types_name)
	print("Dictionaries created")
	data_rel = pd.read_csv(out_rel_name, sep='\t',header = 0)
	data_rel = data_rel.to_numpy(dtype=str)
	return data_rel


def get_reac_prod_dirs(data_dir,file_name): 
	### 
	file_name = data_dir+file_name
	molecules = open(file_name,'r')
	lines = molecules.readlines()
	### building rel_gephi output
	i =0
	reactants = {}
	products = {}
	for line in lines:
		line=line.rstrip('\n')
		a=line.split('\t')

		if a[2] =="1":
			if a[0] in products.keys():
				products[a[0]].append(a[1])
			else:
				products[a[0]]=[a[1]]
		if a[2] =="-1":
			if a[0] in reactants.keys():
				reactants[a[0]].append(a[1])
			else:
				reactants[a[0]]=[a[1]]		
		i += 1

	return reactants,products


def get_spontan_dir(reac_d,equi_data_dir,equilibrator_file_name,is_energy_data):
	equilibrator_file_name = equi_data_dir+equilibrator_file_name
	reactants = reac_d.copy()
	spontan = {}
	if is_energy_data:
		#data = pd.read_csv(equilibrator_file_name, sep='\t',header = 0, dtype = str)
		data = pd.read_csv(equilibrator_file_name, sep='\t',header = None, dtype = str)
		data = data.to_numpy()
		#for Number,rid,name,left,right,dg_kj_mole,dg_error_kj_mole,Error in data: # Set here the format of the enerngy file
		for rid,left,right,dg_kj_mole,dg_error_kj_mole in data: # Set here the format of the enerngy file
			if dg_kj_mole != "NoValue":
				spontan[rid] = float(dg_kj_mole)
			else:
				spontan[rid] = dg_kj_mole		
	else:
		for rea in reactants.keys():
			spontan[rea] = "NoValue"
			#spontan[rea[:-2]] = "NoValue"

	return spontan


def get_reaction_name(data_dir,file_name):
	file_name = data_dir+file_name
	molecules = open(file_name,'r')
	lines = molecules.readlines()
	i =0
	reaction_name_d = {}
	for line in lines:
		line=line.rstrip('\n')
		a=line.split('\t')
		reaction_name_d[a[0]] = a[3]
	return reaction_name_d 



if __name__ == "__main__":

	print("Test 1 :")
	data_dir = "Data/"
	file_name="rels_5.txt"
	react_name_dict = get_reaction_name(data_dir,file_name)
	equi_data_dir = "Equilibrator-Data/"
	
	reaction = "37000_1"
	print("reaction name of ",react_name_dict[reaction])


	print("Test 2 :")
	#equilibrator_file_name = "pH7.4ReactionsDGValuesFinal.txt"  #Put "none" in case we don't have this table yet
	equilibrator_file_name =  "reaction_results.txt"
	is_energy_data = True
	reac_d, prod_d = get_reac_prod_dirs(data_dir,file_name)
	spontaneous_d= get_spontan_dir(reac_d,equi_data_dir,equilibrator_file_name,is_energy_data)
	data_pandas = get_data_pd(data_dir,file_name,spontaneous_d)


	print(data_pandas[0:10])

	reaction = "143938_0"
	#print("DG of the reaction "+ reaction +":",spontaneous_d[reaction[:-2]])
	print("DG of the reaction "+ reaction +":",spontaneous_d[reaction])
	
	print("Test 3 :")
	reaction = "143938_0"
	print("reactants of the reaction "+ reaction +":", reac_d[reaction])
	print("products of the reaction "+ reaction +":", prod_d[reaction])

