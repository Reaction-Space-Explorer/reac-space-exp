from matplotlib import pyplot
import pandas as pd
#import sys
import json
import numpy as np 
import struct
import csv
import ast
from random import sample
import math
import pathlib
from pathlib import Path
import shutil
#sys.setrecursionlimit(20**6) 
from collections import deque
import matplotlib.patches as mpatches
import skunk
from matplotlib.offsetbox import AnnotationBbox
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.offsetbox import AnchoredText
import re 
from collections import Counter

# Use only with the reading_mod_outs library and without the search_autocatalysis 
# and equilibrator_generator libraries.This due to the os use restrictions.


### To do : ########
create_equilibrator_data = True
to_search_autocatalysis = False
to_draw_autocatalysis = False


#### Data parameters
data_dir = "../main/glucose/Neo4j_Imports/rels/"
file_name = "rels_5.txt"
mod_output_molecules = "glucose_degradation_output.csv"


#### Equilibrator parameters
#equilibrator_file_name = "pH7.4ReactionsDGValuesFinal.txt" #Put "none" in case we don't have this and set the follow parameters in order to get the 
equilibrator_file_name = "reaction_results.txt"
equi_data_dir = "Equilibrator-Data/"
is_energy_data = True
to_reset_cache = False
upload_new_compounts_to_cache = False # Use in True only if you want to reset all your databese. Untils this moment i dont how to call a compound if this has no been created corrrectluy.
set_ph_reac_pass = 10
set_ionic_reac_pass = "100 mM"

#### Autocatalysis search parameters
set_results_dir = "AutoCycles-Files-Stubbs/"
#set_results_dir = "Test_results/"
to_seach = "C(C(CC(C(O)=O)=O)O)(O)=O" #Put "all" in case do you want search for all the nodes
cut_energy  = 0.0 # In order to define the spontaneous reaction DG<cut_energy -> spontaneous reaction
erase_water = True

### Drawing parameters 
draws_directory = "AutoCycles-Draws-Stubs/"
only_spontaneous = True



if __name__ == "__main__":

	if create_equilibrator_data:
		import equilibrator_generator_2 as equi
		import reading_mod_outs as re_mod
		reac_dict,prod_dict = re_mod.get_reac_prod_dirs(data_dir,file_name)
		equi.run_equilibrator(equi_data_dir,data_dir,mod_output_molecules,to_reset_cache,upload_new_compounts_to_cache,reac_dict,prod_dict,set_ph_reac_pass,set_ionic_reac_pass)
		print("Wait until the compunts.sqlite is complete created and the compunts.sqlite-journal has disapeared")
		exit()


	if to_search_autocatalysis:
		import search_autocatalysis_3 as search_auto
		import reading_mod_outs as re_mod
		reac_dict,prod_dict = re_mod.get_reac_prod_dirs(data_dir,file_name)
		spontaneous_d= re_mod.get_spontan_dir(reac_dict,equi_data_dir,equilibrator_file_name,is_energy_data).copy()
		data_pandas = re_mod.get_data_pd(data_dir,file_name,spontaneous_d)
		search_auto.search_autocatalysis_cycles (set_results_dir,spontaneous_d,data_pandas,cut_energy, erase_water,to_seach)


	if to_draw_autocatalysis:
		import draw_cycles_14 as draw_cy
		import reading_mod_outs as re_mod
		reac_dict,prod_dict = re_mod.get_reac_prod_dirs(data_dir,file_name)
		spontaneous_d= re_mod.get_spontan_dir(reac_dict,equi_data_dir,equilibrator_file_name,is_energy_data).copy()
		react_name_dict = re_mod.get_reaction_name(data_dir,file_name)
		data_pandas = re_mod.get_data_pd(data_dir,file_name,spontaneous_d)
		draw_cy.draw_cycles(draws_directory,only_spontaneous,reac_dict,prod_dict,spontaneous_d,react_name_dict,set_results_dir)


