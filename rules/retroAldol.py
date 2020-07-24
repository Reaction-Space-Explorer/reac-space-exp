retroAldol = [ruleGMLString("""rule [
	ruleID "Retro Aldol"
	left [
		edge [ source 3 target 4 label "-" ]
		edge [ source 4 target 5 label "-" ]
		edge [ source 5 target 6 label "-" ]
	]
	context [
		node [ id 1 label "C" ]
		node [ id 2 label "O" ]
		node [ id 3 label "C" ]
		node [ id 4 label "C" ]
		node [ id 5 label "O" ]
		node [ id 6 label "H" ]
		edge [ source 1 target 2 label "=" ]
		edge [ source 1 target 3 label "-"]
	]
	right [
		edge [ source 4 target 5 label "=" ]
		edge [ source 3 target 6 label "-" ]
	]
	constrainAdj [ id 4 op "=" count 0 
		 nodeLabels [ label "N" label "S" ]
			edgeLabels [ label "-" label "="]
	]
	constrainAdj [ id 1 op "=" count 0 
		nodeLabels [ label "N" label "S" label "O" ]
		edgeLabels [ label "-" ]
	]
]
""")]