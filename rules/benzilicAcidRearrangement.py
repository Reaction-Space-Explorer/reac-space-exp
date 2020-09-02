benzilicAcidRearrangement = ruleGMLString("""rule [
	ruleID "Benzilic Acid Rearrangement"
	labelType "term"
	left [
		edge [ source 1 target 2 label "-" ]
		edge [ source 3 target 6 label "=" ]
		# break H20 -> OH and H
		edge [ source 8 target 9 label "-" ]
	]
	context [
		node [ id 1 label "*" ]
		node [ id 2 label "C" ]
		node [ id 3 label "C" ]
		node [ id 4 label "*" ]
		node [ id 5 label "O" ]
		node [ id 6 label "O" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 3 target 4 label "-" ]
		edge [ source 2 target 5 label "=" ]

		# and an H2O
		node [ id 7 label "H" ]
		node [ id 8 label "O" ]
		node [ id 9 label "H" ]
		edge [ source 7 target 8 label "-" ]
	]
	right [
		edge [ source 1 target 3 label "-" ]
		edge [ source 3 target 6 label "-" ]
		edge [ source 6 target 9 label "-" ]
		edge [ source 8 target 2 label "-" ]
	]
	# neither carbonyl should be a part of a carboxylic acid, amide, etc.
	constrainAdj [ id 2 op "=" count 0
		nodeLabels [ label "N" label "O" label "S" ]
		edgeLabels [ label "-" ]
	]
	constrainAdj [ id 3 op "=" count 0
		nodeLabels [ label "N" label "O" label "S" ]
		edgeLabels [ label "-" ]
	]
]
""")


benzilicAcidRearrangement_inv = ruleGMLString("""rule [
	ruleID "Benzilic Acid Rearrangement (inverse)"
	labelType "term"
	left [
		edge [ source 1 target 3 label "-" ]
		edge [ source 3 target 6 label "-" ]
		edge [ source 6 target 9 label "-" ]
		edge [ source 8 target 2 label "-" ]
	]
	context [
		node [ id 1 label "*" ]
		node [ id 2 label "C" ]
		node [ id 3 label "C" ]
		node [ id 4 label "*" ]
		node [ id 5 label "O" ]
		node [ id 6 label "O" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 3 target 4 label "-" ]
		edge [ source 2 target 5 label "=" ]

		# and an H2O
		node [ id 7 label "H" ]
		node [ id 8 label "O" ]
		node [ id 9 label "H" ]
		edge [ source 7 target 8 label "-" ]
	]
	right [
		edge [ source 1 target 2 label "-" ]
		edge [ source 3 target 6 label "=" ]
		# break H20 -> OH and H
		edge [ source 8 target 9 label "-" ]
	]
	# should not form an amide or thioic acid
	constrainAdj [ id 2 op "=" count 0
		nodeLabels [ label "N" label "S" ]
		edgeLabels [ label "-" ]
	]
	constrainAdj [ id 3 op "=" count 0
		nodeLabels [ label "N" label "S" ]
		edgeLabels [ label "-" ]
	]
	# carbon 3 should not already be attached to a carboynl, imine etc. group
	constrainAdj [ id 3 op "=" count 0
		nodeLabels [ label "O" label "N" label "S" ]
		edgeLabels [ label "=" ]
	]
]
""")