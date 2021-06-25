# The purpose of this is to provide an "inverse" rule to elimination

# Prefer producing the keto form instead of the enol form (avoid tautomers)
hydration1 = ruleGMLString("""rule [
	ruleID "Hydration of C(=O)C"
	left [
		edge [ source 1 target 2 label "=" ]
		edge [ source 3 target 4 label "-" ]
		edge [ source 5 target 6 label "-" ]
	]
	context [
		node [ id 1 label "C" ]
		node [ id 2 label "O" ]
		node [ id 3 label "C" ]
		node [ id 4 label "H" ]
		#H2O
		node [ id 5 label "O" ]
		node [ id 6 label "H" ]
		node [ id 7 label "H" ]

		edge [ source 1 target 3 label "-" ]
		edge [ source 5 target 7 label "-" ]
	]
	right [
		edge [ source 1 target 2 label "-" ]
		edge [ source 1 target 4 label "-" ]
		edge [ source 2 target 6 label "-" ]
		edge [ source 5 target 3 label "-" ]
	]
	constrainAdj [ id 1 op "=" count 0 
		edgeLabels [ label "-" ]
		nodeLabels [ label "N" label "S" label "O" ]
	]
	constrainAdj [ id 3 op "=" count 0
		edgeLabels [ label "=" label "-" ] # avoid gem diols and forming carboxylic acids, amides
		nodeLabels [ label "O" label "N" label "S" ]
	]
	# make sure it doesn't form enols; for the case of C=C, the suceeding rule is more useful
	constrainAdj [ id 3 op "=" count 0
		nodeLabels [ label "C" ]
		edgeLabels [ label "=" ]
	]
]
""")


hydration2 = ruleGMLString("""rule [
	ruleID "Hydration of C=C(O)"
	left [
		edge [ source 1 target 2 label "=" ]
		edge [ source 3 target 4 label "-" ]
	]
	context [
		node [ id 1 label "C" ]
		node [ id 2 label "C" ]

		node [ id 3 label "H" ]
		node [ id 4 label "O" ]
		node [ id 5 label "H" ]
		edge [ source 4 target 5 label "-" ]
	]
	right [
		edge [ source 1 target 2 label "-" ]
		edge [ source 4 target 2 label "-" ] # join -OH
		edge [ source 3 target 1 label "-" ] # join -H
	]
	constrainAdj [ id 2 op "=" count 0
		nodeLabels [ label "N" label "S" ]
		edgeLabels [ label "-" label "=" ]
	]
	# It should not form gem-diols or aminols
	constrainAdj [ id 2 op "=" count 0
		nodeLabels [ label "O" label "N" label "S" ]
		edgeLabels [ label "-" ]
	]
]
""")

# irreversible
enamine_hydration = ruleGMLString("""rule [
	ruleID "Enamine Hydration and Elimination"
	labelType "term"
	left [
		edge [ source 1 target 2 label "=" ]
		edge [ source 2 target 3 label "-" ]

		edge [ source 7 target 6 label "-" ] # break O-H bonds
		edge [ source 7 target 8 label "-" ]
	]
	context [
		node [ id 1 label "C" ]
		node [ id 2 label "C" ]
		node [ id 10 label "*" ]
		node [ id 11 label "*" ]
		edge [ source 1 target 10 label "-" ]
		edge [ source 1 target 11 label "-" ]

		node [ id 3 label "N" ]
		edge [ source 3 target 4 label "-" ]
		node [ id 4 label "*" ]
		edge [ source 3 target 5 label "-" ]
		node [ id 5 label "*" ]

		node [ id 6 label "H" ]
		node [ id 7 label "O" ]
		node [ id 8 label "H" ]
	]
	right [
		edge [ source 1 target 2 label "-" ] # C-C
		edge [ source 1 target 6 label "-" ] # add H to C
		edge [ source 2 target 7 label "=" ] # C=O
		edge [ source 8 target 3 label "-" ] # H-N
	]
	# Rs Should not be a part of C=O or C=N
	constrainAdj [ id 4 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" ]
	]
	constrainAdj [ id 5 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" ]
	]
	constrainAdj [ id 10 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" ]
	]
	constrainAdj [ id 11 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" ]
	]
]""")
