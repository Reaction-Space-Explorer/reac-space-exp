import matplotlib.pyplot as plt
import math
import re 
import sys
import glob
import os
import os.path
from pathlib import Path
import shutil
import numpy as np
import matplotlib.colors as mcolor
import skunk
from matplotlib.offsetbox import AnnotationBbox
import math
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.offsetbox import AnchoredText

def CountFrequency(my_list):
	freq = {}
	for item in my_list:
		if (item in freq):
			freq[item] += 1
		else:
			freq[item] = 1
	return freq

def plot_histogram(mol,results_dir,plots_dir):
	list_of_files = glob.glob(results_dir+"AutoCatCycles*")

	filed = Path(plots_dir)
	if filed.exists ():
		question = input("Do you want to erase the previos " + plots_dir + " directory (yes/not):")
		if question == "yes" or question == "Yes" or question == "YES":
			shutil.rmtree(plots_dir)
			os.mkdir(plots_dir)
	else:
		os.mkdir(plots_dir)

	for file_name_out in list_of_files:
		len_root = len(results_dir) + len("AutoCatCycles_")
		molecule =file_name_out[len_root:-4] 
		print(file_name_out)
		print(len_root)
		if mol != "all":
			if  molecule != mol:
				continue

		print("Hereeeee")
		file_name = results_dir + "AutoCatCycles_" + molecule + ".txt"

		filed = Path(file_name)
		if not filed.is_file():
			print("There isn't the file: ",file_name)
			exit()	

		len_cycle = []
		len_cycle_spo = []
		cycles = open(file_name,'r')
		lines = cycles.readlines()
		for line in lines:
			line=line.rstrip('\n')
			line=line.split('\t')
			spo = line[-1]
			l0 = int(line[3])/2
			l1 = int(line[5])/2
			l2 = int(line[7])/2
			l = l0+ max(l1,l2)
			len_cycle.append(l)
			if spo == "Strictly_Sponstaneous":
				len_cycle_spo.append(l)

		fig= plt.figure(1, figsize=(8,8))
		fig.clf()
		ax = fig.add_subplot(111)

		labels =CountFrequency(len_cycle).keys()
		values =CountFrequency(len_cycle).values()		
		ax.bar(labels,values,log = True, color ="black",width = 0.4)
		ax.set_xticks(range(1,20))
		ax.set_ylim( (10**-1,10**5) )
		ax.set_xlim( (0,20) )
		ax.set_xlabel("Cycle length",fontsize = 12, fontweight='bold')
		ax.set_ylabel('Frequency',fontsize = 12, fontweight='bold')
		fig.tight_layout()
		svgs = skunk.insert({})
		with open(plots_dir + molecule + "_total" + ".svg", 'w') as fa:
		    fa.write(svgs)
		#plt.savefig(plots_dir + molecule + "_total" +  ".png",transparent=True,dpi=800)
		#plt.show()

		fig= plt.figure(1, figsize=(8,8))
		fig.clf()
		ax = fig.add_subplot(111)
		
		labels =CountFrequency(len_cycle_spo).keys()
		values =CountFrequency(len_cycle_spo).values()		
		ax.bar(labels,values,log = True, color ="black",width = 0.4)
		ax.set_xticks(range(1,20))
		ax.set_ylim( (10**-1,10**5) )
		ax.set_xlim( (0,20) )
		ax.set_xlabel("Cycle length",fontsize = 12, fontweight='bold')
		ax.set_ylabel('Frequency',fontsize = 12, fontweight='bold')
		fig.tight_layout()
		svgs = skunk.insert({})
		with open(plots_dir + molecule + "_spontaneous" + ".svg", 'w') as fa:
		    fa.write(svgs)	
		#plt.savefig(plots_dir + molecule + "_spontaneous" + ".png",transparent=True,dpi=800)
		#plt.show()

		fig= plt.figure(1, figsize=(8,8))
		fig.clf()
		ax = fig.add_subplot(111)

		labels =CountFrequency(len_cycle).keys()
		values =CountFrequency(len_cycle).values()		
		ax.bar(labels,values,log = True, color ="black",width = 0.4)
		ax.set_xticks(range(1,20))
		ax.set_ylim( (10**-1,10**5) )
		ax.set_xlim( (0,20) )
		ax.set_xlabel("Cycle length",fontsize = 12, fontweight='bold')
		ax.set_ylabel('Frequency',fontsize = 12, fontweight='bold')

		box = skunk.Box(250,250,"subplot")
		ab = AnnotationBbox(box,(15,10**3.5),frameon =False,xybox= (0,0),xycoords='data',boxcoords='offset points',box_alignment=(+0.5,+0.5))
		ax.add_artist(ab)
		
		code = "obabel -:"+"'"+molecule+"'"+" -xb none -a -O "+"'"+plots_dir+molecule+".svg"+"'"
		os.system(code)		

		box_2 = skunk.Box(200,200,"subplot_2")
		ab_2 = AnnotationBbox(box_2,(15.5,10**0.9),frameon =False,xybox= (0,0),xycoords='data',boxcoords='offset points',box_alignment=(+0.5,+0.5))
		ax.add_artist(ab_2)
		
		#svgs = skunk.insert({"subplot" : plots_dir + molecule + "_spontaneous" + ".svg" })
		#sku_dir={}
		#sku_dir["subplot"] = plots_dir+molecule+".svg"
		#sku_dir["subplot_2"] = plots_dir + molecule + "_spontaneous" + ".svg"
		text = "A"
		ax.annotate(text, xy=(1, 10**4.5), xytext=(0, 0), fontsize=30,
			    xycoords='data', textcoords='offset points',
			    bbox=dict(facecolor='white', alpha=0.0),
			    horizontalalignment='left', verticalalignment='center',zorder = 1000)
		text = "B"
		ax.annotate(text, xy=(12.8, 10**4.5), xytext=(0, 0), fontsize=30,
			    xycoords='data', textcoords='offset points',
			    bbox=dict(facecolor='white', alpha=0.0),
			    horizontalalignment='left', verticalalignment='center',zorder = 1000)
				
		svg = skunk.insert({ "subplot_2": plots_dir+molecule+".svg", "subplot": plots_dir + molecule + "_spontaneous" + ".svg"})	
		with open(plots_dir + molecule + "_mixed" + ".svg", 'w') as fa:
		    fa.write(svg)
		


def plot_feeder_incycle(dir_feeders,plots_dir,subplot_dir,length,text_to):

	filed = Path(plots_dir+subplot_dir)
	if filed.exists ():
		question = input("Do you want to erase the previos " + plots_dir + " directory (yes/not):")
		if question == "yes" or question == "Yes" or question == "YES":
			shutil.rmtree(plots_dir+subplot_dir)
			os.mkdir(plots_dir+subplot_dir)
	else:
		os.mkdir(plots_dir+subplot_dir)		

	for key in list(dir_feeders.keys())[0:length]:
		print(key)
		code = "obabel -:"+"'"+key+"'"+" -xb none -a -O "+"'"+plots_dir+subplot_dir+key+".svg"+"'"
		os.system(code)
	fig= plt.figure(1, figsize=(14,8))
	fig.clf()
	ax = fig.add_subplot(111)
	ax.bar(list(range(1,length+1)),list(dir_feeders.values())[0:length],log = False, color ="black",width = 0.4)
	
	i = 0
	skunk_dir = {}
	for key in list(dir_feeders.keys())[0:length]:
		i += 1
		box = skunk.Box(100,100,"subplot_"+key)
		ab = AnnotationBbox(box,(i,dir_feeders[key]),frameon =False,xybox= (0,0),xycoords='data',boxcoords='offset points',box_alignment=(+0.5,0))
		ax.add_artist(ab)
		skunk_dir["subplot_"+key] = plots_dir+subplot_dir+key+".svg"
		

	ax.annotate(text_to, xy=(0.85*length, 1.3*list(dir_feeders.values())[0]  ), xytext=(0, 0), fontsize=24,
		    xycoords='data', textcoords='offset points',
		    bbox=dict(facecolor='white', alpha=0.2),
		    horizontalalignment='left', verticalalignment='center',zorder = 1000)
			    
	ax.set_xticks(range(1,length+1))
	ax.set_ylim( (0,1.4*list(dir_feeders.values())[0] ))
	#ax.set_xlim( (0,20) )
	ax.set_xlabel("Molecules",fontsize = 12, fontweight='bold')
	ax.set_ylabel('Frequency',fontsize = 12, fontweight='bold')	
	#plt.xticks(rotation = 90)	
	svg = skunk.insert(skunk_dir)	
	with open(plots_dir + subplot_dir+ text_to + ".svg", 'w') as fa:
		fa.write(svg)

	#plt.show()
		

if __name__ == "__main__":
	molecule = "all"
	set_results_dir = "AutoCycles-Files/"
	plot_dir = "Plots_dir/"
	subplot_dir = "Feeders_dir/"
	figure_caption = "Feeders"
	plot_histogram(molecule,set_results_dir,plot_dir)
	

