import pandas as pd
import os

from pandas.core.frame import DataFrame

"""
A few utility methods to extract useful information from the dump files.

In particular, currently I have put a method that reads the 'rels' files in the Neo4J 
"""

# this folder contains the "reactions" (called 'rels' as they represent 'edges' between nodes (molecules))
rels_dir = 'glucose/Neo4j_Imports/rels/'

class ReactionObject:
	def __init__(self, rxn_id, sub_id, sources, targets):
		self.rxn_id = rxn_id
		self.sub_id = sub_id
		self.sources = sources
		self.targets = targets
	def getNumSources(self):
		num_sources = len(self.sources)
		return numSources
	def getNumTargets(self):
		num_targets = len(self.targets)
		return num_targets

def read_all_gens():
	"""
	Returns: a list sorted by generation (increasing order) containing dataframe objects, each
	containing the reactions of each generation. Note that the rows in each dataframe are unique
	(there is no repetition of G1 and G2 rows in G3 dataframes)
	"""
	# the names of 
	column_names = ['REACTION_ID', 'MOLECULE_SMILES', 'CONSUMPTION/CREATION', 'RULE_NAME']

	# A list of pandas dataframes
	all_dataframes = []
	file_names = sorted(os.listdir(rels_dir))
	for f_name in file_names:
		# the last character of file name (excluding extension) is the generation
		gen = int(f_name.rstrip('.txt')[-1]) 
		rel_data = pd.read_csv(rels_dir + f_name, sep='\t', names=column_names)
		# now because of the way those Neo4J imports were written, gen 3 has data of gen 2 and 1
		# since i'm already reading the files in an "order", I know I need to delete rows in this
		# df that are already in the dataframes that I appended to 'all_dataframes'
		for prev_gen_df in all_dataframes:
			# keep only rows that aren't already there in prev gens
			rel_data = rel_data[~rel_data.isin(prev_gen_df)]
			# this creates NaNs
			rel_data.dropna(inplace=True)
			#rel_data_clean.reset_index(drop=True, inplace=True) # I didn't have a specific reason for this
			# but I thought in case I use indices at some point I would want them reset
		all_dataframes.append(rel_data)
		#print(rel_data)
	return all_dataframes


def get_count_map(dataframe_list):
	"""
	"""
	generation = 1 # assuming the dataframe list starts at G1
	# I shall store the counts in a 'map' {gen_int: {count_dict}} where the 'dict' is a dictionary containg
	# the counts of each rule by generation
	# TODO: replace this giant 'map' with a DataFrame
	# exporting will then just become a one-line thing.
	gen_count_map = {}
	# list of all rules that have been used in the chemistry (only those that have appeared in the
	# dictionary)
	used_rules_list = []
	# now iterate over all dataframes
	for data in dataframe_list:
		# a dictionary of the format {'RULE_NAME': count}
		rule_count_dict = {}
		# Note: since multiple rows in the dataframe comprise a single reaction
		# you need to count all entries of '4_1' as a single rule
		last_id = '0_0'
		for i in range(len(data)):
			# split '4_1' into '4' and '1'
			reac_id = data.iloc[i]['REACTION_ID']
			# if the same reaction id is repeated, this is not a diff reaction
			if reac_id == last_id: 
				pass
			else:
				last_id = reac_id
				rule = data.iloc[i]['RULE_NAME']
				if rule in rule_count_dict.keys():
					rule_count_dict[rule] += 1
				else: # if an entry for that rule isn't there already
					rule_count_dict[rule] = 1
		# add any rules to the list which are not in there already
		for rule in rule_count_dict.keys():
			if rule not in used_rules_list:
				used_rules_list.append(rule)
		# put the whole dictionary in that map
		gen_count_map[generation] = rule_count_dict
		print('Generation:', generation)
		print(rule_count_dict)
		generation += 1 # increment before loop switches to next generation's data
	return gen_count_map, used_rules_list


def find_rule_count(output_file):
	"""
	"""
	dataframe_list = read_all_gens()
	gen_count_map, used_rules_list = get_count_map(dataframe_list)
	# Now create a new dataframe
	table = DataFrame({'Rule': used_rules_list})
	#print(table)
	for gen, count_dict in gen_count_map.items():
		# so in a given generation, not all reactions will have been used and
		# now store count of reactions in the same sequence as that in 'used_rules_list'
		all_rxns_count = []
		for rule in used_rules_list:
			if rule in count_dict.keys():
				all_rxns_count.append(count_dict[rule])
			else:
				all_rxns_count.append(0)
		table[f'Generation {gen}'] = all_rxns_count
	# Sort the table in descending order
	total_gens = len(gen_count_map) # the number of generations
	print(f'Exported reaction count of {total_gens} gens')
	table.sort_values(by=f'Generation {total_gens}', ascending=False, inplace=True)
	table.to_csv(output_file, sep='\t', index=False)


def count_edges_gen():
	"""
	The way the counting is being done here is the following:
	Say, A+B+C -> D+E occurs via 2 different reaction rules. There will be edges from A->D, A->E,
	B->D, B->E, and so on. Each of them is counted as a distinct node.

	The sum of number of edges across generations should turn out to be equal to
	sum([e.numSources * e.numTargets for e in dg.edges]). I have tested this for the glucose degradation
	network and it turns out to work as expected.
	"""
	dataframe_list = read_all_gens()
	edge_count_list = []
	'''edges_per_gen = []
	last_rxn_id = -1
	num_sources = num_targets = 0
	for data in dataframe_list:
		# don't count multiple reactions associated with the same edges twice. Only subscript 0 should be counted
		if '_0' not in data['REACTION_ID']:
			continue
		else:
			if last_rxn_id == :'''
	for data in dataframe_list: # over all generations
		last_rxn_id = ''
		edge_count = 0
		sources = [] # list of source nodes for a given reaction
		targets = [] # list of target nodes for a given reaction
		for i in range(len(data)): # loop over all entries
			entry = data.iloc[i]
			# The reaction ids are going to be strings here. They haven't been converted into ints at this point
			rxn_id, sub_id  = entry['REACTION_ID'].split('_') # i know it looks sloppy but it's not as important
			#print(f'Rxn id: {rxn_id}, sub_id = {sub_id}')
			if rxn_id != last_rxn_id and last_rxn_id != '':
				#print(f'Done with rxn_id: {last_rxn_id}')
				# We have finished loading a reaction from entries, now let's create an object out of it
				#reac = ReactionObject(last_rxn_id, sub_id, sources, targets)
				# NOTE: The following line decides how edges are counted.
				edges_in_rxn = len(sources) * len(targets)
				#print('Edges this rxn:', edges_in_rxn)
				# add this number to the total count of edges for this gen
				edge_count += edges_in_rxn
				# now clean the lists to make space for the next reaction in sequence
				sources = []
				targets = []
				last_rxn_id = rxn_id
			# NOTE: the stuff below needs to be made to happen even if 
			# the normal stuff
			if last_rxn_id == '':
				last_rxn_id = rxn_id
			if sub_id == '0':
				if entry['CONSUMPTION/CREATION'] == -1:
					sources.append(entry['MOLECULE_SMILES']) # reactant
				else: 
					targets.append(entry['MOLECULE_SMILES']) # product
		# Once the loop finishes, the last reaction will be left out. Don't forget to add that.	
		#reac = ReactionObject(rxn_id, sub_id, sources, targets)
		edges_in_rxn = len(sources) * len(targets)
		#print('Edges this rxn:', edges_in_rxn)
		edge_count += edges_in_rxn
		print(f'Edges in this gen : {edge_count}')
		edge_count_list.append(edge_count)
	print('Edge count across all gens: ', sum(edge_count_list))


count_edges_gen()
#find_rule_count('hcn_rule_count_dec2020.tsv')