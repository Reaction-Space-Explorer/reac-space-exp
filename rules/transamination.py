transamination = ruleGMLString("""rule [
	ruleID "Transamination" 
    labelType "term"
	left [
		edge [ source 2 target 3 label "=" ]
		edge [ source 9 target 11 label "-" ]
		edge [ source 11 target 12 label "-" ]
	]   
	context [
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 4 label "-" ]
		edge [ source 4 target 5 label "=" ]
		edge [ source 4 target 6 label "-" ]
		edge [ source 6 target 7 label "-" ]
		edge [ source 11 target 13 label "-" ]
		edge [ source 11 target 14 label "-" ]
		edge [ source 8 target 9 label "-" ]
		edge [ source 9 target 10 label "-" ]
		node [ id 1 label "*" ]
		node [ id 2 label "C" ]
		node [ id 3 label "O" ]
		node [ id 4 label "C" ]
		node [ id 5 label "O" ]
		node [ id 6 label "O" ]
		node [ id 7 label "H" ]
		node [ id 8 label "H" ]
		node [ id 9 label "N" ]
		node [ id 10 label "H" ]
		node [ id 11 label "C" ]
		node [ id 12 label "H" ]
		node [ id 13 label "*" ]
		node [ id 14 label "*" ]
	]   
	right [
		edge [ source 2 target 12 label "-" ]
		edge [ source 2 target 9 label "-" ]
		edge [ source 11 target 3 label "=" ]
	]   
]""")

transamination_inv = ruleGMLString("""rule [
	ruleID "Transamination (inverse)" 
    labelType "term"
	left [
		edge [ source 2 target 12 label "-" ]
		edge [ source 2 target 9 label "-" ]
		edge [ source 11 target 3 label "=" ]
	]     
	context [
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 4 label "-" ]
		edge [ source 4 target 5 label "=" ]
		edge [ source 4 target 6 label "-" ]
		edge [ source 6 target 7 label "-" ]
		edge [ source 11 target 13 label "-" ]
		edge [ source 11 target 14 label "-" ]
		edge [ source 8 target 9 label "-" ]
		edge [ source 9 target 10 label "-" ]
		node [ id 1 label "*" ]
		node [ id 2 label "C" ]
		node [ id 3 label "O" ]
		node [ id 4 label "C" ]
		node [ id 5 label "O" ]
		node [ id 6 label "O" ]
		node [ id 7 label "H" ]
		node [ id 8 label "H" ]
		node [ id 9 label "N" ]
		node [ id 10 label "H" ]
		node [ id 11 label "C" ]
		node [ id 12 label "H" ]
		node [ id 13 label "*" ]
		node [ id 14 label "*" ]
	]   
	right [
		edge [ source 2 target 3 label "=" ]
		edge [ source 9 target 11 label "-" ]
		edge [ source 11 target 12 label "-" ]
	]
	constrainAdj [
		id 11 op "=" count 0
		# should not be a part of a carboxylic acid, amide, etc.
		nodeLabels [ label "O" label "N" label "S" ]
		edgeLabels [ label "-" ]
	]
]""")
