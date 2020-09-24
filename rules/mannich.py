mannich_addition = ruleGMLString("""rule [
	ruleID "Mannich Addition"
	labelType "term"
	left [
		edge [ source 7 target 10 label "-" ]
		edge [ source 12 target 15 label "-" ]
	]
	context [
		node [ id 1 label "*" ]
		node [ id 2 label "C" ]
		edge [ source 1 target 2 label "-" ]
		node [ id 3 label "O" ]
		edge [ source 2 target 3 label "=" ]
		node [ id 4 label "C" ]
		edge [ source 2 target 4 label "-" ]
		node [ id 5 label "H" ]
		edge [ source 4 target 5 label "-" ]
		node [ id 6 label "*" ]
		edge [ source 4 target 6 label "-" ]
		node [ id 7 label "C" ]
		edge [ source 4 target 7 label "-" ]
		node [ id 8 label "*" ]
		edge [ source 7 target 8 label "-" ]
		node [ id 9 label "H" ]
		edge [ source 7 target 9 label "-" ]
		node [ id 10 label "O" ]
		node [ id 11 label "H" ]
		edge [ source 10 target 11 label "-" ]
		node [ id 12 label "N" ]
		node [ id 13 label "*" ]
		edge [ source 12 target 13 label "-" ]
		node [ id 14 label "*" ]
		edge [ source 12 target 14 label "-" ]
		node [ id 15 label "H" ]
	]
	right [
		edge [ source 12 target 7 label "-" ]
		edge [ source 15 target 10 label "-" ]
	]
	constrainAdj [
		id 13 op "=" count 1
		nodeLabels [ label "N" ]
		edgeLabels [ label "-"]
	]
	constrainAdj [
		id 13 op "=" count 0
		nodeLabels [ label "O" ]
		edgeLabels [ label "="]
	]
	constrainAdj [
		id 14 op "=" count 1
		nodeLabels [ label "N" ]
		edgeLabels [ label "-" ]
	]
	constrainAdj [
		id 14 op "=" count 0
		nodeLabels [ label "O" ]
		edgeLabels [ label "=" ]
	]
]
""")

# Can't be inverted automatically due to the constrainAdj so here's redundant code
mannich_addition_rev = ruleGMLString("""rule [
	ruleID "Mannich Addition (reverse)"
	labelType "term"
	left [
		edge [ source 12 target 7 label "-" ]
		edge [ source 15 target 10 label "-" ]
	]
	context [
		node [ id 1 label "*" ]
		node [ id 2 label "C" ]
		edge [ source 1 target 2 label "-" ]
		node [ id 3 label "O" ]
		edge [ source 2 target 3 label "=" ]
		node [ id 4 label "C" ]
		edge [ source 2 target 4 label "-" ]
		node [ id 5 label "H" ]
		edge [ source 4 target 5 label "-" ]
		node [ id 6 label "*" ]
		edge [ source 4 target 6 label "-" ]
		node [ id 7 label "C" ]
		edge [ source 4 target 7 label "-" ]
		node [ id 8 label "*" ]
		edge [ source 7 target 8 label "-" ]
		node [ id 9 label "H" ]
		edge [ source 7 target 9 label "-" ]
		node [ id 10 label "O" ]
		node [ id 11 label "H" ]
		edge [ source 10 target 11 label "-" ]
		node [ id 12 label "N" ]
		node [ id 13 label "*" ]
		edge [ source 12 target 13 label "-" ]
		node [ id 14 label "*" ]
		edge [ source 12 target 14 label "-" ]
		node [ id 15 label "H" ]
	]
	right [
		edge [ source 7 target 10 label "-" ]
		edge [ source 12 target 15 label "-" ]
	]
	constrainAdj [
		id 13 op "=" count 1 # One N would be the amine it itself is attached to
		nodeLabels [ label "N" ]
		edgeLabels [ label "-"]
	]
	constrainAdj [
		id 13 op "=" count 0
		nodeLabels [ label "O" ]
		edgeLabels [ label "="]
	]
	constrainAdj [
		id 14 op "=" count 1
		nodeLabels [ label "N" ]
		edgeLabels [ label "-" ]
	]
	constrainAdj [
		id 14 op "=" count 0
		nodeLabels [ label "O" ]
		edgeLabels [ label "=" ]
	]
]
""")