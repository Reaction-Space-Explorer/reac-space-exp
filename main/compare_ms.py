import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem import MolFromSmiles, MolToSmiles


glucose_data = pd.read_csv(os.path.join("..", 'data/Glucose_MeOH.csv'), delimiter=',',
							 skiprows=[0,1,2,3,4,5,6,7,8, 1039, 1040], usecols=['Mass', 'Kendrick Mass'])
# lines 1039 and 1040 in glucose_data and 1179, 1180 in dextrose_data are "***********"s
# so I had to skip them. First 9 lines contain comments
dextrose_data = pd.read_csv(os.path.join('..', 'data/Dextrose_MeOH.csv'), delimiter=',',
							 skiprows=[0,1,2,3,4,5,6,7,8, 1179, 1180], usecols=['Mass', 'Kendrick Mass'])

# Drop rows containing NaNs
glucose_data.dropna(axis=0, how='any')
dextrose_data.dropna(axis=0, how='any')

# I made a dictionary, whose purpose will become obvious shortly
obs_masses = {}

# The +1.007276 in the mass was done to account for a proton, as this is negative ESI data
for row in glucose_data.itertuples(index=True, name='Pandas'):
	# The CSV was a bit weird, ignore the [0], I just had to access a particular element of the tuple
	obs_masses.update({row[0][1] + 1.007276: 'False'}) # [1] refers to the Mass
	# I set the flag to 'False' to mean that it hasn't been found in our simulation
for row in dextrose_data.itertuples(index=True, name='Pandas'):
	obs_masses.update({getattr(row, 'Mass') + 1.007276: 'False'})

# Mass of simulated species calculated by RDKit might differ a bit from experimental obs
# The mass is that of negative ion, so we need to make a correction for an electron
error_margin = 0.01


def compare_sims(dg, generation, print_extra=False):
	"""
	Compare the fraction of the simulated molecules molecules that are present in MS data

	Keyword arguments:
	smiles_list -- a list of the smiles strings present in the 
	generation -- which generation is this (is it the 2nd or some nth round of glucose hydrolysis?)
	print_extra -- whether or not to print molecules (or reactions forming them) 
				which aren't in the MS data
	"""
	vertices = [v for v in dg.vertices]
	smiles_list = [v.graph.smiles for v in dg.vertices]
	# RDKit Mol object for each smiles string in the network
	simulated_mols = [MolFromSmiles(smiles) for smiles in smiles_list]
	# Calculate exact weight of each molecule
	weights = [ExactMolWt(mol) for mol in simulated_mols]
	# wt: frequency of that weight
	simulated_weights = {}
	# wt: found in data
	simulations_found = {}
	# Mol object: associated weight
	molecule_mass_dict = {}
	# The mass spectra can cover molecules only with mass > 150 amu
	for i in range(len(simulated_mols)):
		if weights[i] >= 150.0:
			# initially, set the value to 'False'. False means it hasn't yet been found in the simulations
			simulations_found[weights[i]] = False
			# calculate the frequency of each weight
			simulated_weights[weights[i]] = weights.count(weights[i])
			molecule_mass_dict.update({simulated_mols[i]: weights[i]}) # .iloc[i,0]
	# Now check which of MOD's generated molecules are in the MS data
	# Update the flag for those observations which have been found to match our simulations
	for obs in obs_masses.keys():
		for sim in simulated_weights.keys():
			if abs(sim - obs) <= error_margin:
				#print(f'Observed mass: {obs} matches simulated structure with weight {sim}')
				obs_masses[obs] = True
				simulations_found[sim] = True
				break
	# no of masses in the data that have been found
	observed_count = 0
	for val in obs_masses.values():
		if val == True:
			observed_count += 1
	print()
	# Elements which are flagged 'False' aren't in the intersection set of simulation and experiment
	matching_sims = 0
	for val in simulations_found.values():
		if val == True:
			matching_sims += 1
	# Now print this info
	print(f'In round {generation}, {observed_count} of {len(obs_masses)} data points (about {100*observed_count/len(obs_masses.keys())}%) in the MS were found in the simulation')
	print(f"The reaction network has {len(weights)} total structures, among these there are {len(simulated_weights.keys())} unique masses")
	print(f'In round {generation}, {matching_sims} of {len(simulated_weights.keys())} ({100*matching_sims/len(simulated_weights.keys())}%) of unique masses in MOD generated structures were seen in MS data')
	print(f'Out of total {len(obs_masses.keys())} points in the MS data, {observed_count} match MOD generated structures')

	# To see which molecules' masses are not in the mass spectra
	if print_extra == True:
		for i in range(len(weights)):
			if weights[i] >= 150 and simulations_found[weights[i]] == False:
				print(f"Molecule {smiles_list[i]} with weight {weights[i]} not in mass spectra")
				postSection("Structure")
				vertices[i].graph.print()
				postSection("Reaction")
				dg2 = DG(graphDatabase=inputGraphs)
				with dg2.build() as b:
					in_edges = [edge for edge in vertices[i].inEdges]
					for edge in in_edges:
						d = Derivations()
						d.left = [source.graph for source in edge.sources]
						d.rules = [rule for rule in edge.rules]
						d.right = [target.graph for target in edge.targets]
						b.addDerivation(d)
					print(f"Printing reaction for non-matching data point: {d}")
				dg2.print()

def make_mass_spectra(smiles_list):
	molecules = [MolFromSmiles(smiles) for smiles in smiles_list]
	weights = [ExactMolWt(mol) for mol in molecules]
	highest_mass = max(weights)
	least_mass = min(weights)
	# make a bar graph of the masses simulated by MOD.
	plt.hist(weights, bins=range(500))
	plt.xlabel("Exact Mass")
	plt.ylabel("Frequency")
	plt.title("Mass spectra of the molecules simulated in the reaction network.")
	plt.show()

if False:
	plt.hist(obs_masses.keys(), bins=range(len(obs_masses.keys())))
	plt.title("The mass spectral data itself")
	plt.xlabel("Exact mass")
	plt.xlim(0, max(obs_masses.keys())+100)
	# TODO: take relative abundance of the molecules into account.
	# plt.ylabel("Frequency")
	plt.show()
	