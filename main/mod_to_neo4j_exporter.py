
def neo4j_export_new(dg_obj, smiles_output_path):
	"""
	A method to create exports of 'relations' (reactions) for Neo4j using the dumps of a network. This involves retrieving
	which reaction happened in which generation. To do that, the smiles txt output is made use of

	NOTE: This method shall include reactions/nodes that existed in say, generations G1 and G2
	 in the output of G3 as well. (This feature is not new and was present in the previous dumps as well).
	"""

	# a dictionary which holds a list of smiles of the species which appeared for a given generation (starting from 1)
	# make sure the smiles output txt was the one that belonged to the same DG dump.
	gen_smiles_map = read_smiles_output(smiles_output_path)

	# the dictionary above will hold the list of smiles that appeared only in a particular generation (say, G3)
	# and not the ones before it (i.e. G1 and G2 smiles won't be in there). I needed the "cumulative" list
	cumulative_species_list = []
	for gen, smiles_list in gen_smiles_map.items():
		# add this generation's smiles to the list
		cumulative_species_list.extend(smiles_list)
		with open(f'Neo4j_Imports/rels/rels_{gen}.txt', 'a') as rels_dump:
			for edge in dg_obj.edges:
				source_smi = [s.graph.smiles for s in edge.sources]
				target_smi = [t.graph.smiles for t in edge.targets]
				# to tell if this reaction corresponds to this particular generation, it should be sufficient to know
				# if all 'target' species first appeared in this gen (or before, b/c say, H2O may have been out there early on)
				happened_this_gen = True
				for target in target_smi:
					if target not in cumulative_species_list:
						happened_this_gen = False
				if happened_this_gen:
					# NOTE: This assumes only rule was associated with this reaction, which may not be true
					for source in source_smi:
						rels_dump.write(source + '\t' + str(edge.id) + '\t' + list(edge.rules)[0].name + '\n')
					for target in target_smi:
						rels_dump.write(str(edge.id) + '\t' + target + '\t' + list(edge.rules)[0].name + '\n')
		
		with open(f'Neo4j_Imports/nodes/nodes_{gen}.txt', 'a') as nodes_dump:
			for v in dg.vertices:
				if v.graph.smiles in cumulative_species_list:
					nodes_dump.write(str(v.graph.id) + '\t' + v.graph.smiles + '\t' + str(v.graph.exactMass) + '\t' +
						 'Molecule' + '\n')


def export_to_neo4j(dg_obj, generation_num):
	"""
	Get export files for importing into Neo4j by iterating across the hypergraph's
	nodes (molecules) and rels/edges (edges).
	"""
	dg_files_path = "Neo4j_Imports"
	generation_num = str(generation_num)
	
	# write nodes txt file
	with open(dg_files_path + '/nodes/nodes_' + generation_num + '.txt', 'a') as nodes_file:
		for v in dg_obj.vertices:
			#print(v)
			# node_id, smiles_str, exact_mass, node_label #Romulo: termolecular species Md and radiation Hf has ficticious mass
			if (v.graph.smiles != "[Hf]" and v.graph.smiles != "[Md]"):
				nodes_file.write(str(v.id) + "," + v.graph.smiles + "," + str(v.graph.exactMass) + ",Molecule" "\n")
			
			if (v.graph.smiles == "[Hf]" ): 
				nodes_file.write(str(v.id) + "," + v.graph.smiles + "," + "Without mass" + ",Radiation" "\n")
			
			if (v.graph.smiles == "[Md]" ): 
				nodes_file.write(str(v.id) + "," + v.graph.smiles + "," + "Mass of the most abundat" + ",Termolecular" "\n")
					
	# write relationships (rels) text file	
	with open(dg_files_path + "/rels/rels_" + generation_num + ".txt", 'a') as rels_file:
		print("#################### Inicio de dg.edges#####################")
		print (str(dg_obj.edges))
		print("#################### Fin de dg.edges#####################")
		for e in dg_obj.edges:
			print("edge id: " + str(e.id))
			print("For neo4j:")
			print(e)
			for r in e.rules: # e.rules has only one element
				print("\tRule:" + str(r)+ "," +str(r.id))
			for s in e.sources:
				print("\tSource:" + str(s) + "," + str(s.graph.smiles)+"," + str(s.id))
			for t in e.targets:
				print("\tTarget:" + str(t) + "," + str(t.graph.smiles)+"," + str(t.id))
			
			
			# only one target per edge, but is collected in a list; same with reaction rule
			#target = list(e.targets)[0]
			rule = list(e.rules)[0]
			#print(str(e.rules))
			
			for target in e.targets:
			# write line for each target molecule
				for source in e.sources:
					source_smiles = str(source.graph.smiles)
					target_smiles = str(target.graph.smiles)
					reaction_rule = str(rule)
					# edge_id, source_smiles_str, target_smiles_str, reaction_rule_str
					rels_file.write(str(e.id) + "," + source_smiles + "," + target_smiles + "," + reaction_rule + "\n")
