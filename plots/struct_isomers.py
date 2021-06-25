"""
Count cumulative structural isomers as a function of generation
Plot N_isomers(gen) vs. exact mass 
"""

from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem import MolFromSmiles

import matplotlib.pyplot as plt
import matplotlib.gridspec as gspec


def count_struct_isomers(smiles_list):
	"""
	Counts the number of molecules with the same molecular formula
	Keyword arguments:
	smiles_list -- a list of smiles strings of the set/subset of molecules to look at
	Returns: 
	"""
	# formula: isomer count
	dict_isomers = {}
	# formula : smiles list
	dict_smiles = {}
	# weight : isomer count
	dict_exactwt = {}
	
	for mol_smiles in smiles_list:
		mol = MolFromSmiles(mol_smiles)
		formula = CalcMolFormula(mol)
		weight = ExactMolWt(mol)
		if formula in dict_isomers.keys():
			dict_isomers[formula] += 1 # increase the isomer count by 1
			dict_smiles[formula].append(mol_smiles) # These are MOD's smiles, not RDKit's
			dict_exactwt[weight] += 1 # Weight calculated by RDKit, not MOD's in-built
		else:
			dict_isomers[formula] = 1
			dict_smiles[formula] = [mol_smiles]
			dict_exactwt[weight] = 1
	return dict_exactwt # modify this as per your needs


# This has now become redundant
"""
def plot_spectra(exactwt_freq_dict, gen):
	'''
	Stacked Bar Plot using matplotlib
	'''
	# which marker (dot) color to use for which generation
	#gen_colors = ['#ffa600', '#ff6361', '#bc5090', '#58508d', '#003f5c']
	
	#
	#gen_colors= ['gray', 'gold', 'teal', 'dodgerblue', 'orangered']
	gen_colors = ["royalblue", "forestgreen", 'gold', 'crimson', 'black'] #"darkorchid"
	
	weights = exactwt_freq_dict.keys()
	freqs = exactwt_freq_dict.values()
	plt.vlines(x=weights, ymin=0, ymax=freqs, color=list(reversed(gen_colors))[gen-1], label=f'Generation {gen}')
	#plt.bar(weights, freqs, width=0.5, color=list(reversed(gen_colors))[gen-1], label=f'Generation {gen}')
	plt.xlabel("Exact Mass")
	plt.ylabel("Cumulative Frequency")
	plt.yscale("log")
	#ax[gen-1].set_yscale('log')
	#ax[gen-1].set_yticks([10,100,1000])
	#plt.xlim(l right=202)
	plt.ylim(bottom=0.625)
	plt.title("Mass spectra of simulated glucose network")
	#plt.setp(stemlines, linestyle="-", color=gen_colors[gen-1], linewidth=0.5)
	#plt.setp(markers, markersize=2, label=f"Generation {gen}")
	#ax[gen-1].legend(loc='upper left')
	#plt.show()"""


def plot_lollipop(exactwt_freq_dict, gen, shared_axis=True):
	"""
	Make a (shared axis) lollipop plot of the spectrum.
	"""
	# which marker (dot) color to use for which generation
	gen_colors = ["blue", "red", "cyan", "green", "magenta"]
	
	weights = list(exactwt_freq_dict.keys())
	# normalized freqs
	freqs = list(exactwt_freq_dict.values())
	# custom ticks (avoid overlapping, increase tick range for certain subplots, etc.)
	axis_ticks = [
		[10, 100],
		[1,10,100],
		[10,100],
		[1,10,100,1000],
		[1,10,100,1000]
	]
	if shared_axis == True:
		# if basefmt is not " " it will draw a coloured horizontal baseline at y=0
		(markers, stemlines, baseline) = axes[gen-1].stem(weights, freqs, basefmt=" ",
				 markerfmt=f"ko", use_line_collection=True) # replace 'k' by {gen_colors[gen-1][0]} for color by geneartion
		#axes[gen-1].bar(weights, freqs, color='black', width=0.25, label=f'Generation {gen}')
		axes[gen-1].set_yscale('log')
		axes[gen-1].set_yticks([10, 100, 1000]) # equal ticks for all subplots
		# use below line for diff ticks for each subplot.
		#axes[gen-1].set_yticks(axis_ticks[gen-1])
		plt.setp(stemlines, linestyle="-", color='gray', linewidth=0.75)
		plt.setp(markers, markersize=1, label=f"Generation {gen}")
		axes[gen-1].legend(loc='upper left')
		# get rid of space between subplots
		plt.subplots_adjust(wspace=0, hspace=0)
	else:
		(markers, stemlines, baseline) = plt.stem(weights, freqs, basefmt="gray",
				 markerfmt=f"{gen_colors[gen-1][0]}o", use_line_collection=True )
		plt.setp(stemlines, linestyle="-", color='gray', linewidth=0.75)
		plt.setp(markers, markersize=2, label=f"Generation {gen}")
		#plt.title("Mass spectra of simulated glucose network")
		plt.yscale("log")
		plt.legend()
	
	

# Only one fig
#fig = plt.figure(figsize=(8,8))
#ax = fig.add_subplot(111)
fig, axes = plt.subplots(5, figsize=(8,8), sharex=True, sharey=True)
#fig.suptitle("Mass spectra of the model glucose reaction network", y=0.91)

# Plot exact wt vs. number of compounds
with open("../main/glucose/glucose_degradation_output_10mar.txt") as output:
	lines = output.readlines()
	gen_smiles_dict = {}
	for line in lines:
		comps = line.split("\t")
		gen = int(comps[0][1])
		## Don't add G0, that's just initial reactants.
		if gen == 0:
			continue
		smiles_str = comps[1]
		if gen in gen_smiles_dict:
			gen_smiles_dict[gen].append(smiles_str)
		else:
			gen_smiles_dict[gen] = [smiles_str]
	
	# List of smiles upto Nth generation
	cumulative_list = []
	for gen, smiles_list in gen_smiles_dict.items():
		cumulative_list.extend(smiles_list)
	count = 0
	# use this for normalization later on
	highest_peak = -1
	for gen in range(len(gen_smiles_dict.keys()), 0, -1):
		exactwt_count = count_struct_isomers(cumulative_list)
		#print(exactwt_count)
		#plot_spectra(exactwt_count, gen)
		#plot_subplots(exactwt_count, gen)
		highest_peak = max(highest_peak, max(exactwt_count.values()))
		print(f'highest freq: {highest_peak}')
		plot_lollipop(exactwt_count, gen)
		#ax2 = fig.add_subplot(111)
		#plot_lollipop(exactwt_count, gen)
		cumulative_list = [x for x in cumulative_list if x not in gen_smiles_dict[gen]]
	#plt.legend(loc='upper left')
	plt.xlabel('Exact Mass')
	plt.ylabel('Cumulative Frequency')
	plt.tight_layout()
	plt.savefig('struct_isomer_freq_gray.jpg', dpi=300)
	plt.show()
