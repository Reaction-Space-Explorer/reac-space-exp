import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FuncFormatter
import datetime
import time

from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

molecules = open('glucose_degradation_output.csv','r')
lines = molecules.readlines()
counter = 0
with open('YSP_ReactionNetwork_Descriptor_Data_November2021.csv', 'w') as the_file:
	the_file.write("Generation,Id,NumStereoIsomers"+'\n') 	
	for line in lines:
		counter +=1
		line=line.rstrip('\n')
		line=line.split('\t')
		m = Chem.MolFromSmiles(line[1])
		isomers = tuple(EnumerateStereoisomers(m))
		numste = str(len(isomers))
		the_file.write(line[0]+","+line[1]+","+numste+'\n') 

molecules = open('YSP_ReactionNetwork_Descriptor_Data_November2021.csv','r')
lines = molecules.readlines()
isom = []
degree = []
gen=[]
for line in lines[1:]:
	line=line.rstrip('\n')
	line=line.split(',')
	gen.append(int(line[0][1:]))
	isom.append(int(line[2]))
	degree.append(int(line[5]))


gen2 = np.array(gen)
isom2 =np.array(isom)
degree2 = np.array(degree)


fig, ax1 = plt.subplots()
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
ax1.set_xscale("log")
ax1.set_xlabel('Degree', fontsize = 10)
ax1.set_ylabel('Stereoisomers',fontsize = 10)
scatter = ax1.scatter(degree2, isom2, c= gen2, label = "Stereoisomers", s=10)


legend1 = ax1.legend(*scatter.legend_elements(),loc="upper right", title="Generation")
ax1.add_artist(legend1)


#h1, l1 = ax1.get_legend_handles_labels()
#ax1.legend(h1, l1, loc='upper left',fontsize = 14,ncol=1,handlelength=0, labelspacing=0.15,markerscale=1,frameon=False)
fig.tight_layout()
fig.savefig('produc.png',dpi=400)