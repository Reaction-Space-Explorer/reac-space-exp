include("common.py")

knoevenagel_c = [ruleGMLString("""rule [
	ruleID "Knoevenagel C"
	labelType "term"
	left [
		edge [ source 1 target 2 label "=" ]
		edge [ source 3 target 4 label "-" ]
	]
	context [
		node [ id 1 label "C" ]
		node [ id 2 label "O" ]

		node [ id 3 label "C" ]
		node [ id 4 label "H" ]
		node [ id 5 label "*" ]
		node [ id 6 label "C" ]
		node [ id 7 label "C" ]

		edge [ source 3 target 6 label "-" ]
		edge [ source 3 target 7 label "-" ]
		edge [ source 3 target 5 label "-" ]
	]
	right [
		edge [ source 1 target 3 label "-" ]
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 4 label "-" ]
	]
	# The C=O that merges should not be a part of a -(C=O)NH2, etc.
	constrainAdj [ id 1 op "=" count 0
		nodeLabels [ label "O" label "N" label "S" ]
		edgeLabels [ label "-" ]
	]
	# The R can be either C#N or C=O
	constrainAdj [ id 6 op "=" count 1
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" label "#" ] # Warning: this may allow =N to be matched as well
	]
	constrainAdj [ id 7 op ">=" count 1
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "-" label "=" label "#" ] # Warning: this may allow =N to be matched as well
	]
]""")]

knoevenagel_h = [ruleGMLString("""rule [
	ruleID "Knoevenagel H"
	labelType "term"
	left [
		edge [ source 1 target 2 label "=" ]
		edge [ source 3 target 4 label "-" ]
	]
	context [
		node [ id 1 label "C" ]
		node [ id 2 label "O" ]

		node [ id 3 label "C" ]
		node [ id 4 label "H" ]
		node [ id 5 label "*" ]
		node [ id 6 label "C" ]

		edge [ source 3 target 6 label "-" ]
		edge [ source 3 target 5 label "-" ]
	]
	right [
		edge [ source 1 target 3 label "-" ]
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 4 label "-" ]
	]
	# The C=O that merges should not be a part of a -(C=O)NH2, etc.
	constrainAdj [ id 1 op "=" count 0
		nodeLabels [ label "O" label "N" label "S" ]
		edgeLabels [ label "-" ]
	]
	# The R can be either C#N or C=O
	constrainAdj [ id 6 op "=" count 1
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" label "#" ] # Warning: this may allow =N to be matched as well
	]
	constrainAdj [ id 3 op ">=" count 1
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "-" label "=" label "#" ] # Warning: this may allow =N to be matched as well
	]
]""")]


'''def knoevenagelGen():
	def attach_CN_CdOA(r, root, offset):
		r.context.extend([
			"# CN or C(=O)A",
			'edge [ source %d target %d label "-" ]' % (root, offset),
			'node [ id %d label "C" ]' % offset,
		])
		rCopy = r.clone()
		rCopy.name += ", CN"
		rCopy.context.extend([
			'node [ id %d label "N" ]' % (offset + 1),
			'edge [ source %d target %d label "#" ]' % (offset, offset + 1),
		])
		yield rCopy
		rCopy = r.clone()
		rCopy.name += ", C(=O)A"
		rCopy.label = "term"
		rCopy.context.extend([
			'node [ id %d label "O" ]' % (offset + 1),
			'edge [ source %d target %d label "=" ]' % (offset, offset + 1),
			'node [ id %d label "*" ]' % (offset + 2),
			'edge [ source %d target %d label "-" ]' % (offset, offset + 2),
		])
		yield rCopy


	r = RuleGen("Knoevenagel")
	r.left.extend([
		'# C-EWG',
		'edge [ source 0 target 1 label "-" ]',
		'# C(=O)',
		'edge [ source 10 target 11 label "=" ]',
	])
	r.context.extend([
		'# C-EWG',
		'node [ id 0 label "C" ]',
		'node [ id 1 label "H" ]',
		'# C=O',
		'node [ id 10 label "C" ]',
		'node [ id 11 label "O" ]',
	])
	r.right.extend([
		'edge [ source 0 target 10 label "-" ]',
		'edge [ source 1 target 11 label "-" ]',
		'edge [ source 10 target 11 label "-" ]',
	])
	# C-EWG
	for r in attach_H_C(r, 0, 100):
		for r in attach_CN_CdOA(r, 0, 200):
			for r in attach_EWG(r, 0, 300):
				# C=O
				r.name += " |"
				for r in attach_2_H_C(r, 10, 400):
					yield r.loadRule()

knoevenagel = [a for a in knoevenagelGen()]'''