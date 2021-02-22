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
				nodes_file.write(str(v.id) + '\t' + v.graph.smiles + '\t' + str(v.graph.exactMass) + '\t' +
						 'Molecule' + '\n')
			if (v.graph.smiles == "[Hf]" ): 
				nodes_file.write(str(v.id) + "\t" + v.graph.smiles + "\t" + "Without mass" + "\tRadiation" "\n")
			
			if (v.graph.smiles == "[Md]" ): 
				nodes_file.write(str(v.id) + "\t" + v.graph.smiles + "\t" + "Mass of the most abundat" + "\t" + "Termolecular" "\n")
					
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
			# for all rules in e.rules
			for i in range(len(e.rules)):
				# use 'i' as a "sub-index"
				for source in e.sources:
					rels_file.write(f'{e.id}_{i}' + '\t' + source.graph.smiles + '\t' + '-1' + '\t' + list(e.rules)[i].name + '\n')
				for target in e.targets:
					rels_file.write(f'{e.id}_{i}' + '\t' + target.graph.smiles + '\t' + '1' + '\t' + list(e.rules)[i].name + '\n')
		