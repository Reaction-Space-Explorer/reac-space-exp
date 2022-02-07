import pandas as pd
import matplotlib.pyplot as plt
import time

from rdkit.Chem import MolFromSmiles
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt

"""
Count structural isomers, molecular formulae and plot these
There's also an option to dump everything to a .tsv

And yes, I made a new script for this (cleaner code)
"""

formose_dat = pd.read_csv('../main/formose/formose_output_may2021.txt', sep='\t',
						names=['Generation', 'SMILES'])

def compute_quantities(dat):
	"""
	Computes the molecular formula, exact molecular weight and adds them as a column to the passed dataframe
	
	dat --- the DataFrame object containing the table of smiles, etc.

	"""
	weights = []
	formulae = []
	print('Starting computation...')
	start_time = time.time()
	for i in range(len(dat)):
		row = dat.iloc[i]
		mol = MolFromSmiles(row['SMILES'])
		weights.append(ExactMolWt(mol))
		formulae.append(CalcMolFormula(mol))
	end_time = time.time()
	print(f'Finished in: {end_time-start_time:0.2f}s')
	dat['Exact Wt'] = weights
	dat['Formula'] = formulae



def plot_isomer_freq():
	"""
	"""
	fig, axes = plt.subplots(6, figsize=(6,6), sharex=True, sharey=True)
	# G0 would be '#9e0142'
	gen_colors = ['#e95c47', '#fdbf6f', '#ffffbe', '#bfe5a0', '#54aead', '#5e4fa2']
	# To make things cumulative, I'll keep adding prev gens stuff
	prev_gens_data = []

	for gen, gen_data in formose_dat.groupby('Generation'):
		if gen == 'G0': # ignore G0
			continue
		gen_num = int(gen[1]) # 'G1' --> 1

		cumulative_data = pd.concat([*prev_gens_data, gen_data]) # concat to get cumulative
		prev_gens_data.append(gen_data) # add this gen's data to the list
		counts_df = cumulative_data.groupby('Exact Wt').count()

		# Ignore the fact that the y-axis in the following line is 'Formula'. Apparently, all other columns
		# will contain the frequency after .count() is called
		(markers, stemlines, baseline) = axes[gen_num-1].stem(counts_df.index, counts_df['Formula'], basefmt=" ",
				 markerfmt=f"ko", use_line_collection=True)
		plt.setp(stemlines, linestyle="-", color='gray', linewidth=0.75)
		plt.setp(markers, markersize=3.5, mfc= gen_colors[gen_num-1], mec='k', mew=0.4, label=f"Generation {gen_num}")
		#axes[gen_num-1].scatter(counts_df.index, counts_df['Formula'], label=f'Generation {gen_num}',
		#			c=gen_colors[gen_num-1])

		axes[gen_num-1].set_yscale('log')
		axes[gen_num-1].set_yticks([10, 100, 1000])

		axes[gen_num-1].minorticks_on()
		axes[gen_num-1].tick_params(axis='both', labelsize=13)
		axes[gen_num-1].legend(fontsize=13, markerscale=2.2, handletextpad=0.1, loc='upper left')
		#print(counts_df.head())
	plt.xlabel('Exact Mass', fontsize=13)
	axes[3].set_ylabel('Cumulative Frequency', fontsize=13)
	plt.tight_layout()
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.savefig('struct_isomer_count.png', dpi=300)
	plt.show()


compute_quantities(formose_dat)
print(formose_dat.head())
plot_isomer_freq()