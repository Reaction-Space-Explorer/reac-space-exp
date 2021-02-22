"""
I needed to dump the list of reactions that happened in a certain format, so that it could be used
to calculate thermodynamic properties. This file contains the code which will read the txt dumps which
contain the information about the reactions
"""

import os

folder_root = '../main/glucose/Neo4j_Imports'

# this folder contains files that just have the list of nodes (molecules)
nodes_files = os.listdir(folder_root + '/nodes')
# this one contains "relations" between nodes (?reactions?)
rels_files = os.listdir(folder_root + '/rels')

# Creating a container for a reaction makes it less messy to do the processing.

class Reaction:
	"""
	A tiny object to hold a reaction, with info gathered from Jessica's dumps
	"""


	def __init__(self):
		self.rxn_id = None
		self.reactants = []
		self.products = []
		self.name = None


	def __repr__(self):
		"""
		The string representation of a Reaction object. Helpful in printing the reaction in a readable fashion.
		The representation should be

		r_1 + r_2 + .. r_n -> p_1 + p_2 + .. p_m <tabspace> (<reaction_name>)

		where r_i are the reactant smiles and p_j are the product smiles
		"""
		str_out = ' + '.join(self.reactants)
		str_out += ' -> '
		str_out += ' + '.join(self.products)
		str_out += '\t' + f'({self.name})'
		return str_out


# put all the Reaction objects in a list.
reactions_list = []

# since the files contain cumulative list of reactions, get the most recent item in the list
# it will contain all reactions
with open(folder_root + '/rels/' + sorted(rels_files, reverse=False)[1]) as f:
	lines = f.readlines()
	# load each reaction into a Reaction object
	reaction = Reaction()
	for line in lines:
		line = line.rstrip('\n')
		# contents of the line are separated by ','
		components = line.split(',')
		reac_id = components[0]
		reac_smi = components[1]
		prod_smi = components[2]
		# since the separator is ',' and some reaction names contain that (e.g. 'Hemiacetal Formation for 5 membered rings, inverse')
		# the names themselves may get split
		reac_name = ','.join(components[3:]) # join them back with ','s
		# as it goes line after line, if a different reaction id is encountered, it means the last few lines were
		# all the reaction was.
		if reac_id != reaction.rxn_id:
			# so, add the last reaction to the list (except if it's an empty object)
			if reaction.rxn_id is not None:
				reactions_list.append(reaction)
				print(reaction)
			# now create a new, container for the next reaction
			reaction = Reaction()
			reaction.rxn_id = reac_id
			reaction.name = reac_name
		# add the reactant and product smiles to the list
		reaction.reactants.append(reac_smi)
		reaction.products.append(prod_smi)
	# Due to how the above loop was constructed, the last reaction entry might get left out without
	# getting added. Add it now
	if reaction.rxn_id != None and reaction not in reactions_list:
		reactions_list.append(reaction)
		print(reaction)
	print('Finished reading', f.name)