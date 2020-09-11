knoevenagel_c = ruleGMLString("""rule [
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
]""")

knoevenagel_h = ruleGMLString("""rule [
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
]""")

knoevenagel_c_inv = ruleGMLString("""rule [
	ruleID "Knoevenagel C (inverse)"
	labelType "term"
	left [
		edge [ source 1 target 3 label "-" ]
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 4 label "-" ]
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
		edge [ source 1 target 2 label "=" ]
		edge [ source 3 target 4 label "-" ]
	]
	# The -OH should not be a part of a carboxylic acid
	constrainAdj [ id 1 op "=" count 0
		nodeLabels [ label "O" ]
		edgeLabels [ label "=" ]
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
]""")

knoevenagel_h_inv = ruleGMLString("""rule [
	ruleID "Knoevenagel H (inv)"
	labelType "term"
	left [
		edge [ source 1 target 3 label "-" ]
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 4 label "-" ]
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
		edge [ source 1 target 2 label "=" ]
		edge [ source 3 target 4 label "-" ]
	]
	# The -OH should not be a part of a carboxylic acid
	constrainAdj [ id 1 op "=" count 0
		nodeLabels [ label "O" ]
		edgeLabels [ label "=" ]
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
]""")