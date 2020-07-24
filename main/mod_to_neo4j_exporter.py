

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
