##########################################################################
################ This program is to plot cycles from MOD outputs #########
##########################################################################
##########################################################################

######################### Requirements ###################################
# You need to install skunk 
# pip install skunk
# https://github.com/whitead/skunk
#########################################################################
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import skunk
from matplotlib.offsetbox import AnnotationBbox
import math
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.offsetbox import AnchoredText
import re 
import sys
import glob
import os
import os.path
import pathlib
import shutil
from collections import Counter
import subprocess
import reading_mod_outs as re_mod
import ast
#########################################################################

	
#########################################################################
def mod(a, b):
	ret = a%b;
	if ret>=0:
		return ret
	else:
		return ret + b

def pol2cart(rho, phi): # rho is radius, and phi is angle
    x = rho * math.cos(math.radians(phi))
    y = rho * math.sin(math.radians(phi))
    return [x, y]

def get_principal_poligone(x,y,l,n,branch_p,branch_0,branch_1,branch_2,l_feed,angle_feed,reac_dict_1,prod_dict_1,l_text, angle_text):
	#reac_dict_1 = reac_dict.copy()
	#prod_dict_1 = prod_dict.copy()
	branch_principal = branch_p.copy()
	list_vertices = []
	list_reac_feed_coor = []
	list_reac_feed_mol = []
	list_prod_feed_coor = []
	list_prod_feed_mol = []
	water_number_reac = []
	water_number_prod = []
	text_position_list = []
	for i in range(n):
		angle_node = 90-i*360/n
		x1 = x + pol2cart(l,angle_node)[0]
		y1 = y + pol2cart(l,angle_node)[1]
		list_vertices.append((x1,y1))
		#print(branch_principal)
		if branch_principal[i][:-2].isnumeric():
			reac_out, prod_out,water_count_reac,water_count_prod = get_mol_out_cycle(True,i,branch_principal[i],branch_0,branch_1,branch_2,reac_dict_1,prod_dict_1)
			u,v = get_feder_mol_poligone(x1,y1,angle_node,l_feed,angle_feed)
			w = get_text_position(x1 ,y1, angle_node, l_text, angle_text)
			list_reac_feed_coor.append(u)
			list_reac_feed_mol.append(reac_out.copy())
			list_prod_feed_coor.append(v)
			list_prod_feed_mol.append(prod_out.copy())
			water_number_reac.append(water_count_reac)
			water_number_prod.append(water_count_prod)
			text_position_list.append(w)		
		else:
			list_reac_feed_coor.append((0,0))
			list_reac_feed_mol.append([])
			list_prod_feed_coor.append((0,0))
			list_prod_feed_mol.append([])
			water_number_reac.append(0)
			water_number_prod.append(0)
			text_position_list.append((0,0))				
	return list_vertices,list_reac_feed_coor,list_reac_feed_mol,list_prod_feed_coor,list_prod_feed_mol,water_number_reac,water_number_prod,text_position_list

def get_secundary_poligone(center,angle,r_c,d_s,m,branch_s,branch_0,branch_1,branch_2,l_feed,angle_feed,reac_dict,prod_dict,l_text,angle_text):


	reac_dict_2 = reac_dict.copy()
	prod_dict_2 = prod_dict.copy()

	branch_secundary = branch_s.copy()
	(a,b) = center
	list_vertices = []
	list_reac_feed_coor = []
	list_reac_feed_mol = []
	list_prod_feed_coor = []
	list_prod_feed_mol = []
	water_number_reac = []
	water_number_prod = []
	text_position_list = []
	a = a - d_s*math.sin(math.radians(angle/2))
	b = b - d_s*math.cos(math.radians(angle/2))
	r_s = math.sqrt(pow(d_s+r_c*math.cos(math.radians(angle/2)),2) + pow(r_c*math.sin(math.radians(angle/2)),2) )
	angle_new = 2*(180/math.pi)*math.acos( (-d_s - r_c*math.cos(math.radians(angle/2)))  /r_s)
	for j in range(m):
		theta = -90-angle/2+angle_new/2
		theta_2 = theta - angle_new
		angle_node = theta-j*(theta-theta_2)/(m-1)
		x1 = a + pol2cart(r_s,angle_node)[0]
		y1 = b + pol2cart(r_s,angle_node)[1]
		list_vertices.append((x1,y1))

		if branch_secundary[j][:-2].isnumeric():
			reac_out, prod_out,water_count_reac,water_count_prod = get_mol_out_cycle(False,j,branch_secundary[j],branch_0,branch_1,branch_2,reac_dict_2,prod_dict_2)
			u,v = get_feder_mol_poligone(x1,y1,angle_node,l_feed,angle_feed)
			w = get_text_position(x1 ,y1, angle_node, l_text, angle_text)
			list_reac_feed_coor.append(u)
			list_reac_feed_mol.append(reac_out.copy())
			list_prod_feed_coor.append(v)
			list_prod_feed_mol.append(prod_out.copy())
			water_number_reac.append(water_count_reac)
			water_number_prod.append(water_count_prod)
			text_position_list.append(w)

		else:
			list_reac_feed_coor.append((0,0))
			list_reac_feed_mol.append([])
			list_prod_feed_coor.append((0,0))
			list_prod_feed_mol.append([])
			water_number_reac.append(0)
			water_number_prod.append(0)
			text_position_list.append((0,0))		
	sec_cent = [a,b]

	return list_vertices,sec_cent,list_reac_feed_coor,list_reac_feed_mol,list_prod_feed_coor,list_prod_feed_mol,water_number_reac,water_number_prod,text_position_list

def get_mol_out_cycle(in_principal,position,reac_node,branch_0,branch_1,branch_2,reac_dict_3,prod_dict_3):
	#reac_dict_3 = reac_dict.copy()
	#prod_dict_3 = prod_dict.copy()
	cat_reaction = branch_0.split(",")[-1]
	cat_node = branch_0.split(",")[0]
	if branch_1 !="":
		b_p = branch_0+","+branch_1+","+cat_node 
	else:
		b_p = branch_0+","+cat_node
	branch_principal_1 = b_p.split(",")
	if branch_2 != "":
		b_s = cat_reaction+","+branch_2+","+cat_node
	else:
		b_s = cat_reaction+","+cat_node
	branch_secundary_1 = b_s.split(",")

	reac_list = reac_dict_3[reac_node].copy()
	prod_list = prod_dict_3[reac_node].copy()
	len_b0 = len(branch_0.split(","))
	if position != (len_b0-1) and in_principal :
		segment_list_reac = [branch_principal_1[position-1]]
		segment_list_prod = [branch_principal_1[position+1]]
	if position != 0 and not in_principal:
		segment_list_reac = [branch_secundary_1[position-1]]
		segment_list_prod = [branch_secundary_1[position+1]] 
	if (position == (len_b0-1) and in_principal):
		segment_list_reac = [branch_principal_1[position-1]]
		segment_list_prod = [branch_principal_1[position+1],branch_secundary_1.copy()[1]]

	if position == 0 and not in_principal :
		segment_list_reac = [branch_0.split(",")[-2]]
		if branch_1 !="":
			segment_list_prod = [branch_1.split(",")[0],branch_secundary_1.copy()[1]]
		else:
			segment_list_prod = [branch_0.split(",")[0],branch_secundary_1.copy()[1]]	
	for seg in segment_list_reac:
		if seg in reac_list:
			reac_list.remove(seg)
	for seg in segment_list_prod:
		if seg in prod_list:
			prod_list.remove(seg)
	reac_out = reac_list.copy()
	prod_out = prod_list.copy()
	water_count_reac = reac_out.count("O")
	water_count_prod = prod_out.count("O")	
	reac_out = [x for x in reac_out if x != "O"]
	prod_out = [x for x in prod_out if x != "O"]

	return reac_out,prod_out,water_count_reac,water_count_prod

def get_feder_mol_poligone(x,y,angle,l_feed,angle_feed): # n and m are the size of the primary and secundary cicles.
	ang_reac = angle + angle_feed
	ang_prod = angle - angle_feed
	x1_reac = x + pol2cart(l_feed,ang_reac)[0]
	y1_reac = y + pol2cart(l_feed,ang_reac)[1]
	x1_prod = x + pol2cart(l_feed,ang_prod)[0]
	y1_prod = y + pol2cart(l_feed,ang_prod)[1]
	return (x1_reac,y1_reac),(x1_prod,y1_prod)

def get_text_position(x,y,angle,l_text,angle_tex): # n and m are the size of the primary and secundary cicles.
	ang_position = angle + angle_tex
	x1_text = x + pol2cart(l_text,ang_position)[0]
	y1_text = y + pol2cart(l_text,ang_position)[1]
	return (x1_text,y1_text)
	
	
def draw_molecules(branch_l):
	branch_list = branch_l.copy()
	for node in branch_list:
		if not node[:-2].isnumeric():
			code = "obabel -:"+"'"+node+"'"+" -xb none -a -O "+"'"+node+'.svg'+"'"
			#print(code)
			os.system(code)				
def erase_molecules():
	list_of_files = glob.glob('*.svg') 
	for file_name in list_of_files:
		try:
			f = open(file_name)
			os.remove(file_name)
		except IOError:
		    print("File not accessible")
	
###############################################################################
def draw_cycles(draws_directory,only_spontaneous,reac_dict,prod_dict,spontaneous_d,react_name_dict,dir_cycle_results): 
	### Check if we want to erase previous draws
	reac_dict_1 = reac_dict.copy()
	prod_dict_1 = prod_dict.copy()
	spontaneous_d_1 = spontaneous_d.copy()
	react_name_dict_1 = react_name_dict.copy()
	draws_directory = draws_directory
	os.getcwd()
	filed = pathlib.Path(draws_directory)
	if filed.exists ():
		question = input("Do you want to erase the previos AutoCycles-Draws directory (yes/not):")
		if question == "yes" or question == "Yes" or question == "YES":
			shutil.rmtree(draws_directory)
			os.mkdir(draws_directory)
		else:
			print("Finishing the script")
			exit()	
	else:
		os.mkdir(draws_directory)
		

	#######################################
	######### Parameters ##################
	principal_nodes_size = 9.0
	reaction_nodes_size = 1.0
	draw_cycles = True
	figure_size = 70
	factor_princ_len_to_radious = 1.9
	factor_secundary_cycle_radiuos  = 2.7

	node_mole_color = "gold"
	transparence_mole_node = 0.4
	node_reaction_color = "darkgoldenrod"
	transparence_reaction_node = 1 

	arrow_color = "dimgray"
	arrow_width = 0.5
	arrow_headwidth = 4.0
	arrow_headlength = 4.0


	frameon_set = False
	mol_size_porcentage = 50


	l_text = -3
	angle_text = 0
	
	len_feeders = 14
	angle_feeders = 35
	feeder_arrow_color = "lightgray"
	feeder_node_color = "magenta"
	transparence_feeder_node = 0.2 
	##############################################################################
	##############################################################################
	list_of_files = glob.glob(dir_cycle_results+"AutoCatCycles*")
	#list_of_files = glob.glob('Test_results/AutoCatCycles*')
	#print(list_of_files)
	
	for file_name in list_of_files:
		count = 0
		cata_name =file_name[len(dir_cycle_results):-4]
		#print(cata_name)
		os.mkdir("AutoCycles-Draws/"+cata_name)
		#print("AutoCycles-Draws/"+cata_name)
		######################################################################
		with open(file_name) as f:
			lines = f.readlines()
			for line in lines:
				if not line.startswith("#"):
					count += 1
					###############################################################################
					line=line.rstrip('\n')
					line = line.split("\t")
					cat_node = line[0]
					cat_reaction = line [1]
					branch_0 = line[2]
					len_0 = int(line[3])
					branch_1 = line[4]
					len_1 =int(line[5])
					branch_2 = line[6]
					len_2 = int(line[7])
					total_energy = line[8]
					reaction_list = line [9]
					reaction_energy_list = line [10]


					############
					if only_spontaneous and (line[11] == "Non_Strictly_Sponstaneous" or line[11] == "NonValue"):
						continue
					############
					if branch_1 !="":
						b_p = branch_0+","+branch_1
					else:
						b_p = branch_0

					branch_principal = b_p.split(",")
					if branch_2 != "":
						b_s = cat_reaction+","+branch_2+","+cat_node
					else:
						b_s = cat_reaction+","+cat_node
					branch_secundary = b_s.split(",")
					############################################################################
	
					############################################################################
					principal_len_cycle = len_0+len_1
					angle_cat_reaction = 360*(len_0-1)/principal_len_cycle
					secundary_len_cycle = len_2+2
					############################################
					dir_mol_princ = {}
					dir_mol_sec = {}
					idx =0
					for nod in branch_principal:
						dir_mol_princ[idx] = nod
						idx += 1
					idx =0
					for nod in branch_secundary:
						dir_mol_sec[idx] = nod
						idx += 1			
					#print(branch_principal)
					#print(branch_secundary)
					draw_molecules(branch_principal)
					draw_molecules(branch_secundary)

					print("AutoCycles-Draws/"+cata_name)			
					print(count)
					############################################################################
					radious_cycle = principal_len_cycle*factor_princ_len_to_radious
					##############################################################################
					###################### Draw Autoactalytic Cycles
					###############################################################################			
					if draw_cycles == False:
						exit()			
					###############################################################################
					fig= plt.figure(1, figsize=(8,8))
					fig.clf()
					ax = fig.add_subplot(111)
					##############################################################################
					skunk_dir = {}

					##############################################################################
					list_ver,list_reac_feed_coor,list_reac_feed_mol,list_prod_feed_coor,list_prod_feed_mol,water_number_reac,water_number_prod,text_position_list = get_principal_poligone(0,0,radious_cycle,principal_len_cycle,branch_principal,branch_0,branch_1,branch_2,len_feeders,angle_feeders,reac_dict_1,prod_dict_1,l_text, angle_text)
					indexe = 0
					#print(text_position_list)
					for tup in list_ver:								
						idx = list_ver.index(tup)
						if (idx % 2) == 0:
							el= mpatches.Ellipse(tup, principal_nodes_size, principal_nodes_size, angle=30, alpha=transparence_mole_node,color =node_mole_color)
							ax.add_artist(el)
						else:
							el= mpatches.Rectangle((tup[0]-reaction_nodes_size/2,tup[1]-reaction_nodes_size/2), reaction_nodes_size, reaction_nodes_size, angle=0, alpha=transparence_reaction_node,color =node_reaction_color)						
							ax.add_artist(el)
							# Water molecules
							#print(dir_mol_princ)
							reac_smiles = dir_mol_princ.copy()[indexe]
							#print(reac_smiles)
							#energy_reaction = spontaneous_d_1.copy()[reac_smiles[0:-2]]
							energy_reaction = spontaneous_d_1.copy()[reac_smiles]
							if energy_reaction == "NoValue":
								energy_reaction = "NoValue"
							else:
								energy_reaction = "{:.2f}".format(energy_reaction)
							#print(energy_reaction)						  		
							(c,d) = text_position_list[idx]
							#print(text_position_list)
							#ax.text(c, d, "+"+ str(water_number_reac[idx]) + " H2O \n" + "-"+str(water_number_prod[idx]) + " H2O", color='black', size = 6 , ha='center', va = 'center' )
							if water_number_reac[idx] > 0 and water_number_prod[idx] > 0:
								text = "r:" + reac_smiles + "\n "+ str(energy_reaction)+ "kj.mol$^{-1}$" + "\n" + "-" + str(water_number_reac[idx]) + " H$_2$O \n" + "+"+str(water_number_prod[idx]) + " H2O"
							if water_number_reac[idx] > 0 and water_number_prod[idx] == 0:
								text = "r:" +  reac_smiles + "\n "+ str(energy_reaction)+ "kj.mol$^{-1}$" + "\n" + "-" + str(water_number_reac[idx]) + " H$_2$O"
							if water_number_reac[idx] == 0 and water_number_prod[idx] > 0:
								text = "r:" +  reac_smiles + "\n "+ str(energy_reaction)+ "kj.mol$^{-1}$" + "\n" + "+" +str(water_number_prod[idx]) + " H$_2$O"
							if water_number_reac[idx] == 0 and water_number_prod[idx] == 0:
								text = "r:" +  reac_smiles + "\n "+ str(energy_reaction)+ "kj.mol$^{-1}$"
							ax.annotate(text, xy=(c, d), xytext=(0, 0), fontsize=3.8,
							            xycoords='data', textcoords='offset points',
							            bbox=dict(facecolor='white', alpha=0.0),
							            horizontalalignment='center', verticalalignment='center',zorder = 1000)

						if ((idx-1) % 2) == 0:
							el2 = mpatches.Ellipse(list_ver[idx-1], principal_nodes_size, principal_nodes_size, angle=30, alpha=transparence_mole_node,color = node_mole_color)
							ax.add_artist(el2)
						else:
							el2 = mpatches.Rectangle((list_ver[idx-1][0]-reaction_nodes_size/2,list_ver[idx-1][1]-reaction_nodes_size/2), reaction_nodes_size, reaction_nodes_size, angle=0, alpha=transparence_reaction_node,color = node_reaction_color)
							ax.add_artist(el2)
							
						ax.annotate("",
							    xy=tup, xycoords='data',
							    xytext=list_ver[idx-1], textcoords='data',
							    arrowprops=dict(color=arrow_color,
									    width=arrow_width ,
									    headlength = arrow_headwidth,
									    headwidth = arrow_headlength,
									    patchA=el2,
									    patchB=el,
									    shrinkB=0,
									    shrinkA=0,
									    connectionstyle="arc3,rad=-0.3",
									    ),
							    )




						#########################################################################
						if (idx % 2) != 0:
							reac_node = dir_mol_princ.copy()[indexe]
							num_reactants = len(list_reac_feed_mol[idx])
							num_products = len(list_prod_feed_mol[idx])
							tup_rf = list_reac_feed_coor[idx]
							tup_pf = list_prod_feed_coor[idx]

							el_reac_feed= mpatches.Ellipse(tup_rf, principal_nodes_size, principal_nodes_size, angle=30, alpha=transparence_feeder_node,color =feeder_node_color)					
							if num_reactants > 0:
								ax.add_artist(el_reac_feed)							
							el_prod_feed= mpatches.Ellipse(tup_pf, principal_nodes_size, principal_nodes_size, angle=30, alpha=transparence_feeder_node,color =feeder_node_color)					
							if num_products > 0:
								ax.add_artist(el_prod_feed)			

							if num_reactants > 0:
								ax.annotate("",
									    xy=list_ver[idx], xycoords='data',
									    xytext=tup_rf, textcoords='data',
									    arrowprops=dict(color=feeder_arrow_color,
											    width=arrow_width ,
											    headlength = arrow_headwidth,
											    headwidth = arrow_headlength,
											    patchA=el_reac_feed,
											    patchB=el,
											    shrinkB=0,
											    shrinkA=0,
											    connectionstyle="arc3,rad=-0.3",
											    ),
									    )

							if num_products > 0:
								ax.annotate("",
									    xy=tup_pf, xycoords='data',
									    xytext=list_ver[idx], textcoords='data',
									    arrowprops=dict(color=feeder_arrow_color,
											    width=arrow_width ,
											    headlength = arrow_headwidth,
											    headwidth = arrow_headlength,
											    patchA=el,
											    patchB=el_prod_feed,
											    shrinkB=0,
											    shrinkA=0,
											    connectionstyle="arc3,rad=-0.3",
											    ),
									    )
							draw_molecules(list_reac_feed_mol[idx])
							
							for mol in list_reac_feed_mol[idx]:
								#print(mol)
								#print(list_reac_feed_mol[idx][0])	
								box = skunk.Box(mol_size_porcentage,mol_size_porcentage,mol+"_"+str(indexe)+"_feed_reac")
								ab = AnnotationBbox(box,tup_rf,frameon =frameon_set,xybox= (0,0),xycoords='data',boxcoords='offset points',box_alignment=(+0.5,+0.5))
								#ax.scatter(tup_rf[0],tup_rf[1])
								ax.add_artist(ab)
								#print(dir_mol_princ[indexe]+".svg")
								skunk_dir[mol+"_"+str(indexe)+"_feed_reac"]=mol+".svg"
							draw_molecules(list_prod_feed_mol[idx])
							
							
							for mol in list_prod_feed_mol[idx]:
								#print(mol)
								#print(list_reac_feed_mol[idx][0])	
								box = skunk.Box(mol_size_porcentage,mol_size_porcentage,mol+"_"+str(indexe)+"_feed_prod")
								ab = AnnotationBbox(box,tup_pf,frameon =frameon_set,xybox= (0,0),xycoords='data',boxcoords='offset points',box_alignment=(+0.5,+0.5))
								#ax.scatter(tup_pf[0],tup_pf[1])
								ax.add_artist(ab)
								#print(dir_mol_princ[indexe]+".svg")
								skunk_dir[mol+"_"+str(indexe)+"_feed_prod"]=mol+".svg"
							
						##########################################################################
						if (indexe % 2) == 0:
							#print(dir_mol_princ[indexe]+"_"+str(indexe))
							box = skunk.Box(mol_size_porcentage,mol_size_porcentage,dir_mol_princ[indexe]+"_"+str(indexe))
							ab = AnnotationBbox(box,tup,frameon =frameon_set,xybox= (0,0),xycoords='data',boxcoords='offset points',box_alignment=(+0.5,+0.5))
							ax.add_artist(ab)
							#print(dir_mol_princ[indexe]+".svg")
							skunk_dir[dir_mol_princ[indexe]+"_"+str(indexe)]=dir_mol_princ[indexe]+".svg"			
						indexe += 1
						###########################################################################
					indexe = 1
					list_ver_2,sec_center,list_reac_feed_coor,list_reac_feed_mol,list_prod_feed_coor,list_prod_feed_mol,water_number_reac,water_number_prod,text_position_list = get_secundary_poligone((0,0),angle_cat_reaction,radious_cycle,secundary_len_cycle*factor_secundary_cycle_radiuos ,secundary_len_cycle,branch_secundary,branch_0,branch_1,branch_2,len_feeders,angle_feeders,reac_dict_1,prod_dict_1,l_text, angle_text)
					for tup in list_ver_2[1:]:				
						idx = list_ver_2.index(tup)
						if (idx % 2) != 0:
							el= mpatches.Ellipse(tup, principal_nodes_size, principal_nodes_size, angle=30, alpha=transparence_mole_node,color =node_mole_color)
							ax.add_artist(el)
						if (idx % 2) == 0:
							el= mpatches.Rectangle((tup[0]-reaction_nodes_size/2,tup[1]-reaction_nodes_size/2), reaction_nodes_size, reaction_nodes_size, angle=0, alpha=transparence_reaction_node,color =node_reaction_color)
							ax.add_artist(el)		

							reac_smiles = dir_mol_sec.copy()[indexe]
							#print(reac_smiles)
							#energy_reaction = spontaneous_d_1.copy()[reac_smiles[0:-2]]
							energy_reaction = spontaneous_d_1.copy()[reac_smiles]
							if energy_reaction == "NoValue":
								energy_reaction = "NoValue"
							else:
								energy_reaction = "{:.2f}".format(energy_reaction)
							#print(energy_reaction)						  		
							(c,d) = text_position_list[idx]
							#print(text_position_list)
							#ax.text(c, d, "+"+ str(water_number_reac[idx]) + " H2O \n" + "-"+str(water_number_prod[idx]) + " H2O", color='black', size = 6 , ha='center', va = 'center' )
							if water_number_reac[idx] > 0 and water_number_prod[idx] > 0:
								text = "r:" + reac_smiles + "\n "+ str(energy_reaction)+ "kj.mol$^{-1}$" + "\n" + "-" + str(water_number_reac[idx]) + " H$_2$O \n" + "+"+str(water_number_prod[idx]) + " H2O"
							if water_number_reac[idx] > 0 and water_number_prod[idx] == 0:
								text = "r:" +  reac_smiles + "\n "+ str(energy_reaction)+ "kj.mol$^{-1}$" + "\n" + "-" + str(water_number_reac[idx]) + " H$_2$O"
							if water_number_reac[idx] == 0 and water_number_prod[idx] > 0:
								text = "r:" +  reac_smiles + "\n "+ str(energy_reaction)+ "kj.mol$^{-1}$" + "\n" + "+" +str(water_number_prod[idx]) + " H$_2$O"
							if water_number_reac[idx] == 0 and water_number_prod[idx] == 0:
								text = "r:" +  reac_smiles + "\n "+ str(energy_reaction)+ "kj.mol$^{-1}$"
							ax.annotate(text, xy=(c, d), xytext=(0, 0), fontsize=4.5,
							            xycoords='data', textcoords='offset points',
							            bbox=dict(facecolor='white', alpha=0.0),
							            horizontalalignment='center', verticalalignment='center',zorder = 1000)



						if ((idx-1) % 2) != 0:	
							el2 = mpatches.Ellipse(list_ver_2[idx-1], principal_nodes_size, principal_nodes_size, angle=30, alpha=transparence_mole_node,color = node_mole_color)
							ax.add_artist(el2)
						else:
							el2 = mpatches.Rectangle((list_ver_2[idx-1][0]-reaction_nodes_size/2,list_ver_2[idx-1][1]-reaction_nodes_size/2), reaction_nodes_size, reaction_nodes_size, angle=0, alpha=transparence_reaction_node,color =node_reaction_color)
							ax.add_artist(el2)
									
						ax.annotate("",
							    xy=tup, xycoords='data',
							    xytext=list_ver_2[idx-1], textcoords='data',
							    arrowprops=dict(color=arrow_color,
									    width=arrow_width ,
									    headlength = arrow_headwidth,
									    headwidth = arrow_headlength,
									    patchA=el2,
									    patchB=el,
									    shrinkB=0,
									    shrinkA=0,
									    connectionstyle="arc3,rad=-0.5",
									    ),
							    )


						##########################################################################
						#########################################################################
						if (idx % 2) == 0:
							reac_node = dir_mol_sec.copy()[indexe]
							num_reactants = len(list_reac_feed_mol[idx])
							num_products = len(list_prod_feed_mol[idx])
							tup_rf = list_reac_feed_coor[idx]
							tup_pf = list_prod_feed_coor[idx]
							el_reac_feed= mpatches.Ellipse(tup_rf, principal_nodes_size, principal_nodes_size, angle=30, alpha=transparence_feeder_node,color =feeder_node_color)					
							if num_reactants > 0:
								ax.add_artist(el_reac_feed)
							el_prod_feed= mpatches.Ellipse(tup_pf, principal_nodes_size, principal_nodes_size, angle=30, alpha=transparence_feeder_node,color =feeder_node_color)					
							if num_products > 0:
								ax.add_artist(el_prod_feed)

							if num_reactants > 0:
								ax.annotate("",
									    xy=list_ver_2[idx], xycoords='data',
									    xytext=tup_rf, textcoords='data',
									    arrowprops=dict(color=feeder_arrow_color,
											    width=arrow_width ,
											    headlength = arrow_headwidth,
											    headwidth = arrow_headlength,
											    patchA=el_reac_feed,
											    patchB=el,
											    shrinkB=0,
											    shrinkA=0,
											    connectionstyle="arc3,rad=-0.3",
											    ),
									    )


							if num_products > 0:
								ax.annotate("",
									    xy=tup_pf, xycoords='data',
									    xytext=list_ver_2[idx], textcoords='data',
									    arrowprops=dict(color=feeder_arrow_color,
											    width=arrow_width ,
											    headlength = arrow_headwidth,
											    headwidth = arrow_headlength,
											    patchA=el,
											    patchB=el_prod_feed,
											    shrinkB=0,
											    shrinkA=0,
											    connectionstyle="arc3,rad=-0.3",
											    ),
									    )
							draw_molecules(list_reac_feed_mol[idx])
						
							for mol in list_reac_feed_mol[idx]:
								box = skunk.Box(mol_size_porcentage,mol_size_porcentage,mol+"_"+str(indexe)+"_feed_reac")
								ab = AnnotationBbox(box,tup_rf,frameon =frameon_set,xybox= (0,0),xycoords='data',boxcoords='offset points',box_alignment=(+0.5,+0.5))
								ax.add_artist(ab)
								skunk_dir[mol+"_"+str(indexe)+"_feed_reac"]=mol+".svg"
							draw_molecules(list_prod_feed_mol[idx])
							
							
							for mol in list_prod_feed_mol[idx]:
								box = skunk.Box(mol_size_porcentage,mol_size_porcentage,mol+"_"+str(indexe)+"_feed_prod")
								ab = AnnotationBbox(box,tup_pf,frameon =frameon_set,xybox= (0,0),xycoords='data',boxcoords='offset points',box_alignment=(+0.5,+0.5))
								ax.add_artist(ab)
								skunk_dir[mol+"_"+str(indexe)+"_feed_prod"]=mol+".svg"
						##########################################################################
						indexe = idx
						if (idx % 2) != 0 and idx != list_ver_2.index(list_ver_2[-1]):					
							box = skunk.Box(mol_size_porcentage,mol_size_porcentage,dir_mol_sec[indexe]+"_"+str(indexe))
							ab = AnnotationBbox(box,tup,frameon =frameon_set,xybox= (0,0),xycoords='data',boxcoords='offset points',box_alignment=(+0.5,+0.5))
							ax.add_artist(ab)
							skunk_dir[dir_mol_sec[indexe]+"_"+str(indexe)]=dir_mol_sec[indexe]+".svg"
						indexe += 1
				        ####################################################################################
					x_lim_min = -figure_size+sec_center[0]
					x_lim_max = figure_size+sec_center[0]
					y_lim_min = -figure_size+sec_center[1]
					y_lim_max = figure_size+sec_center[1]
					size_plot = x_lim_max  - x_lim_min					
					
					text = ""
					reac_list_tolist = ast.literal_eval(reaction_list)
					for reaction_id in reac_list_tolist:
						text = text + reaction_id.ljust(11," ") + "\n" 
					
					x_position_legend = x_lim_min * 99.9 / 100
					y_position_legend = y_lim_min * 75 / 100
					

					ax.annotate(text, xy=(x_position_legend, y_position_legend), xytext=(0, 0), fontsize=6.5,
						    xycoords='data', textcoords='offset points',
						    bbox=dict(facecolor='white', alpha=0.0),
						    horizontalalignment='left', verticalalignment='center',zorder = 1000)					

					text = ""
					for reaction_id in reac_list_tolist:
						text = text + ":   " + react_name_dict_1[reaction_id] + "\n" 
					ax.annotate(text, xy=(x_position_legend + size_plot/15 , y_position_legend), xytext=(0, 0), fontsize=6.5,
						    xycoords='data', textcoords='offset points',
						    bbox=dict(facecolor='white', alpha=0.0),
						    horizontalalignment='left', verticalalignment='center',zorder = 1000)					 					
					############################################################################################
					ax.set_xlim(x_lim_min,x_lim_max)
					ax.set_ylim(y_lim_min,y_lim_max)
					plt.tick_params(left = False, right = False , labelleft = False , labelbottom = False, bottom = False)
					fig.tight_layout()
					ax.axis('off')
					###########################################################################################
					svgs = skunk.insert(skunk_dir)
					with open(draws_directory+cata_name+"/"+cata_name+"_"+str(count)+".svg", 'w') as fa:
					    fa.write(svgs)
					plt.close()
					erase_molecules()



if __name__ == "__main__":

	print("Test 1: ")
	draws_directory = "AutoCycles-Draws/"
	only_spontaneous = True
	file_name ="rels_5.txt"
	data_dir = "Data/"
	equi_data_dir = "Equilibrator-Data/"
	set_results_dir = "AutoCycles-Files/"
	#set_results_dir = "Test_results/"
	#equilibrator_file_name = "pH7.4ReactionsDGValuesFinal.txt"  #Put "none" in case we don't have this table yet
	equilibrator_file_name =  "reaction_results.txt"
	is_energy_data = True

	reac_dict,prod_dict = re_mod.get_reac_prod_dirs(data_dir,file_name)
	spontaneous_d= re_mod.get_spontan_dir(reac_dict,equi_data_dir,equilibrator_file_name,is_energy_data).copy()
	react_name_dict = re_mod.get_reaction_name(data_dir,file_name)
	draw_cycles(draws_directory,only_spontaneous,reac_dict,prod_dict,spontaneous_d,react_name_dict,set_results_dir)

