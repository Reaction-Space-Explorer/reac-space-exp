with_formaldehyde=False

# Note: with_formaldehyde is a variable that must be initialized before you include main.py
# or all.py or anything in your job file. Otherwise it will crash.
if with_formaldehyde == True:
	cannizarro_hcho_ox = [ruleGMLString("""rule [
		ruleID "Cannizarro 2, HCHO (oxidation)"
		labelType "term"
		left [
			edge [ source 8 target 7 label "-" ]
			edge [ source 2 target 3 label "=" ]
			edge [ source 6 target 10 label "-" ]
		]   
		context [
			#edge [ source 1 target 2 label "-" ]
			edge [ source 2 target 4 label "-" ]
			#edge [ source 5 target 6 label "-" ]
			edge [ source 6 target 11 label "=" ]
			edge [ source 8 target 9 label "-" ]
			#node [ id 1 label "*" ]
			node [ id 2 label "C" ]
			node [ id 3 label "O" ]
			node [ id 4 label "H" ]
			#node [ id 5 label "*" ]
			node [ id 6 label "C" ]
			node [ id 7 label "H" ]
			node [ id 8 label "O" ]
			node [ id 9 label "H" ]
			node [ id 10 label "H" ]
			node [ id 11 label "O" ]
			# Here, I am going to force one of them to be HCHO and not any other aldehyde
			# to reduce number of reactions.
			node [ id 12 label "H" ]
			edge [ source 12 target 6 label "-" ]
		]	
		right [
			edge [ source 6 target 8 label "-" ]
			edge [ source 2 target 3 label "-" ]
			edge [ source 3 target 10 label "-" ]
			edge [ source 2 target 7 label "-" ]
		]
		constrainAdj[ id 2 op "=" count 0
			nodeLabels [ label "O" label "N" label "S" ]
			edgeLabels [ label "-" ]
		]
	]""")]

	cannizarro_hcho_red = [ruleGMLString("""rule[
		ruleID "Cannizarro 2, HCHO (reduction)"
		labelType "term"
		left [
			edge [ source 8 target 7 label "-" ]
			edge [ source 2 target 3 label "=" ]
			edge [ source 6 target 10 label "-" ]
		]   
		context [
			#edge [ source 1 target 2 label "-" ]
			edge [ source 2 target 4 label "-" ]
			#edge [ source 5 target 6 label "-" ]
			edge [ source 6 target 11 label "=" ]
			edge [ source 8 target 9 label "-" ]
			#node [ id 1 label "*" ]
			node [ id 2 label "C" ]
			node [ id 3 label "O" ]
			node [ id 4 label "H" ]
			#node [ id 5 label "*" ]
			node [ id 6 label "C" ]
			node [ id 7 label "H" ]
			node [ id 8 label "O" ]
			node [ id 9 label "H" ]
			node [ id 10 label "H" ]
			node [ id 11 label "O" ]
			# Here, I am going to force one of them to be HCHO and not any other aldehyde
			# to reduce number of reactions.
			node [ id 12 label "H" ]
			edge [ source 12 target 2 label "-" ]
		]	
		right [
			edge [ source 6 target 8 label "-" ]
			edge [ source 2 target 3 label "-" ]
			edge [ source 3 target 10 label "-" ]
			edge [ source 2 target 7 label "-" ]
		]
		constrainAdj[ id 6 op "=" count 0
			nodeLabels [ label "O" label "N" label "S" ]
			edgeLabels [ label "-" ]
		]
	]""")]

else: 
	def gen_cann_glucose():
		"""
		Generate 2 rules for Cannizarro with glucose as one of the oxidant/reductants
		In one case glucose is oxidized and in the other it is reduced.
		"""
		glucose = smiles("O=CC(O)C(O)C(O)C(O)C(O)", name="glucose", add=False)
		glucose_gml = glucose.getGMLString().split("\n")
		glucose_gml = glucose_gml[1:-2] # get rid of the 'graph [' and last ']'
		for i in range(len(glucose_gml)):
			glucose_gml[i] = glucose_gml[i].replace('\t', '')
		comm_context = [
			'# some aldehyde group that needs reducing',
			'node [ id 101 label "C" ]',
			'node [ id 102 label "O" ]',
			'node [ id 103 label "H" ]',
			'#water',
			'node [ id 104 label "O" ]',
			'node [ id 105 label "H" ]',
			'node [ id 106 label "H" ]',
			'edge [ source 104 target 105 label "-" ]',
			]
		ox_rule = RuleGen("Cannizarro 2, glucose (oxidation)")
		red_rule = RuleGen("Cannizarro 2, glucose (reduction)")

		ox_rule.context.extend(comm_context)
		red_rule.context.extend(comm_context)
		ox_break = [
			'edge [ source 12 target 1 label "-" ]',
			'edge [ source 101 target 102 label "=" ]',
			'edge [ source 104 target 106 label "-" ]',
		]

		ox_rule.left.extend(ox_break)

		ox_rule.right.extend([
			'# reduce the aldehyde',
			'edge [ source 12 target 102 label "-" ]',
			'edge [ source 101 target 102 label "-" ]',
			'edge [ source 101 target 106 label "-" ]',
			'#oxidize glucose by adding -OH',
			'edge [ source 1 target 104 label "-" ]',
		])
		# add the glucose part to the context. Note: no line should in 'left' should be in 'context'
		ox_context_glu = []
		for line in glucose_gml:
			if line not in ox_break:
				ox_context_glu.append(line)
		
		ox_rule.context.extend(ox_context_glu)
		ox_rule.context.extend(['edge [ source 101 target 103 label "-" ]'])

		red_rule_left = [
			'edge [ source 1 target 0 label "=" ]',
			'edge [ source 104 target 106 label "-" ]',
			'edge [ source 103 target 101 label "-" ]'
			]
		red_rule.left.extend(red_rule_left)
		red_context_glu = []
		for line in glucose_gml:
			if line not in red_rule_left:
				red_context_glu.append(line)
		red_rule.context.extend(red_context_glu)
		red_rule.context.extend(['edge [ source 102 target 101 label "=" ]'])
		red_rule.right.extend([
			'edge [ source 103 target 0 label "-"] ',
			'edge [ source 1 target 0 label "-" ]',
			'edge [ source 104 target 101 label "-" ]',
			'edge [ source 1 target 106 label "-" ]'
		])
		constraint = [
			'constrainAdj [ id 101 op "=" count 0 ',
			'\tnodeLabels [ label "O" label "N" label "S" ]',
			'\tedgeLabels [ label "-" ]',
			']'
		]
		ox_rule.constraints.extend(constraint)
		red_rule.constraints.extend(constraint)

		return [ox_rule.loadRule(), red_rule.loadRule()]

	cannizarro_glucose_rule = gen_cann_glucose()

	cannizarro_glucose_rule[0].print()