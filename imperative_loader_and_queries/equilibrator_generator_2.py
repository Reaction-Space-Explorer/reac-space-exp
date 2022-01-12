import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd 
import numpy as np 
import csv
import os
from equilibrator_assets.local_compound_cache import LocalCompoundCache
import reading_mod_outs as re_mod
from pathlib import Path

def run_equilibrator(equilibrator_data_directorio,data_dir,mod_output_molecules,to_reset_cache,upload_new_compounts_to_cache,reactant_dir,prods_dir,set_ph_reac_pass,set_ionic_reac_pass):
	to_reset_cache = to_reset_cache	
	mole_output_file = mod_output_molecules
	data_dir = data_dir
	equi_data_dir = equilibrator_data_directorio
	upload_new_compounts_to_cache = upload_new_compounts_to_cache
	if to_reset_cache:
		lc = LocalCompoundCache()
		lc.generate_local_cache_from_default_zenodo(data_dir+"compounds.sqlite")
		lc.load_cache(equi_data_dir+"compounds.sqlite")
		print("defaul zenodo data base uploaded. Please re run the program setting to_reset_cache= False")	
	else:
		lc = LocalCompoundCache()
		lc.load_cache(equi_data_dir+"compounds.sqlite")
	set_ph_reac = set_ph_reac_pass
	set_ionic_reac = set_ionic_reac_pass # format : "100 mM"

	###########################################################
	########### Read reactions ################################
	if upload_new_compounts_to_cache:
		print("uploading")
		molecules = open(data_dir+mole_output_file,'r')
		lines = molecules.readlines()
		list_moleculas = []
		dir_mol = {}
		count = 0
		for line in lines[:500]:
			count += 1
			line=line.rstrip('\n')
			line=line.split('\t')
			dir_mol[line[1]]=str(count)
			list_moleculas.append([line[1],str(count),str(count)])

		##########################################################

		compound_df = pd.DataFrame(data=list_moleculas,columns=["struct","coco_id", "name"])
		lc.add_compounds(compound_df, mol_format="smiles",bypass_chemaxon=False,save_empty_compounds=True,error_log=equi_data_dir+"compound_creation_log.tsv")
		print("Completed upload of new compounds. Please edit any error in compound_creation_log.csv (ssave as *.csv) and re run the program setting upload_new_compounts_to_cache = False ")
		exit()
	co = 0
	mole_database = open(equi_data_dir+"compound_creation_log.csv",'r')
	lines = mole_database.readlines()
	dir_correc_mol = {}
	for line in lines[1:]:
		co += 1
		print(co)
		line=line.rstrip('\n')
		line=line.split('\t')
		if len(line)<4:
			print("correct this line in the compound_creation_log.tsv: ", str(co))
			exit()
		dir_correc_mol[line[1]]= line[4]
	
	
	print("In order to run correctly \n your compound_creation_log.csv file must contain \n all the molecules in your reactions")
	checking = input("Your compound_creation_log.tsv is completed? (yes, no): " )
	if checking == "no":
		print("finishing the program")
	else:
		print("continuing")
	###########################################################
	###########################################################

	to_erase = Path(equi_data_dir+"reaction_results.txt")
	if to_erase.is_file():
		question = input("Do you want to errase the previous output "+ equi_data_dir+"reaction_results.txt" +" (yes,no)")
		if question == "yes" or question == "Yes" or question == "YES":
			os.remove(to_erase)
		else:
			print("Finishing the script")
			exit()	

	to_erase_2 = Path(equi_data_dir+"to_upload.txt")
	if to_erase_2.is_file() :
		question = input("Do you want to errase the previous output "+ equi_data_dir+"to_upload.txt" +" (yes,no)")
		if question == "yes" or question == "Yes" or question == "YES":
			os.remove(to_erase_2)
		else:
			print("Finishing the script")
			exit()			

	###########################################################
	#This method uses the equilibrator_api and the LocalCompoundCache to enable custom-compound use.

	from equilibrator_api import ComponentContribution, Q_
	from equilibrator_api import Reaction

	cc = ComponentContribution(ccache = lc.ccache)

	###########################################################
	############# Reaction ####################################
	###########################################################
	reacs = reactant_dir
	prods = prods_dir
	to_upload_mols = []	
	with open(equi_data_dir+"reaction_results.txt","w") as file_reactions:				
		reaction_count = 0
		for rea in list(reacs.keys())[0:10]:
			reaction_count += 1
			print("#################")
			print("\t \t "+ rea + "\t \t ")
			print("#################")
			print(str(reacs[rea])+"---"+str(prods[rea]))
			dir_reaction= {}
			count = 1
			reac_count = {i:reacs[rea].count(i) for i in reacs[rea]}
			fail_count = 0
			for reactante in reacs[rea]:
				count +=1
				reactante_compound =lc.get_compounds(reactante,return_fails=True)
				print(reactante_compound)
				if not reactante in dir_correc_mol.keys():
					fail_count+=1
					to_upload_mols.append(reactante)
				
				else:
					if dir_correc_mol[reactante] == "failed":
						print("############ reactivo no encontrado ###############")
						fail_count+=1

				dir_reaction[reactante_compound]=-reac_count[reactante]
			count = 1
			prods_count = {i:prods[rea].count(i) for i in prods[rea]}

			for producto in prods[rea]:
				count+=1
				producto_compound = lc.get_compounds(producto,return_fails=True)
				print(producto_compound)
				if not producto in dir_correc_mol.keys():
					fail_count+=1
					to_upload_mols.append(producto)
				else:
					if dir_correc_mol[producto] =="failed":
						print("############ producto no encontrado ###############")	
						fail_count +=1
				dir_reaction[producto_compound]=prods_count[producto]
			
			if fail_count == 0:
				reaction_group = Reaction(dir_reaction)
				
				standard_dg_prime = cc.standard_dg_prime(reaction_group)
				if not reaction_group.is_balanced():
					print('%s is not balanced' %  reaction_group )	
				cc.p_h = Q_(set_ph_reac)  # set pH
				cc.ionic_strength = Q_(set_ionic_reac_pass)  # set I
				dG_prime = cc.dg_prime(reaction_group)
				dG_prime_value_in_kj_per_mol = dG_prime.value.m_as("kJ/mol")
				dG_prime_error_in_kj_per_mol = dG_prime.error.m_as("kJ/mol")
				file_reactions.write(str(rea)+"\t"+str(reacs[rea])+"\t"+str(prods[rea])+"\t"+str(dG_prime_value_in_kj_per_mol)+"\t"+str(dG_prime_error_in_kj_per_mol)+"\n")				
				
				print(dG_prime_value_in_kj_per_mol)
				print(dG_prime_error_in_kj_per_mol)
			else:
				file_reactions.write(str(rea)+"\t"+str(reacs[rea])+"\t"+str(prods[rea])+"\t"+"NoValue"+"\t"+"NoValue"+"\n")


	##################################
	### Reporting some issues ########
	print(to_upload_mols)
	to_upload_mols = list(set(to_upload_mols))
	if len(to_upload_mols)>0:
		print("There are some molecules that isn't in database: See the file to_upload.txt ")
		with open(equi_data_dir+"to_upload.txt","w") as to_upload:
			for mol in to_upload_mols:
				to_upload.write(mol+"\n")
		

if __name__ == "__main__":

	print("Test Equilibrato 1 :")
	mod_output_molecules = "glucose_degradation_output.csv"
	file_name ="rels_5.txt"
	data_dir = "Data/"
	equi_data_dir = "Equilibrator-Data/"
	reac_dict,prod_dict = re_mod.get_reac_prod_dirs(data_dir,file_name)

	to_reset_cache = False
	upload_new_compounts_to_cache = False  # Use in True only if you want to reset all your databese. Untils this moment i dont how to call a compound if this has no been created corrrectluy.
	set_ph_reac_pass = 10
	set_ionic_reac_pass = "100 mM"


	run_equilibrator(equi_data_dir,data_dir,mod_output_molecules,to_reset_cache,upload_new_compounts_to_cache,reac_dict,prod_dict,set_ph_reac_pass,set_ionic_reac_pass)
		
		
