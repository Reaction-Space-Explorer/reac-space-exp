ketoEnolisation = [ruleGMLString("""rule [
	ruleID "Keto-enol migration twice" # to avoid enols (and tautomers in general)
	left[
		edge [ source 1 target 2 label "=" ]

		edge [ source 3 target 4 label "-" ]
		edge [ source 3 target 6 label "-" ]
	]
	context [
		node [ id 1 label "C" ]
		node [ id 2 label "O" ]

		node [ id 3 label "C" ]
		node [ id 4 label "O" ]
		node [ id 5 label "H" ]
		node [ id 6 label "H" ]
		edge [ source 4 target 5 label "-" ]
		edge [ source 1 target 3 label "-" ]
	]
	right [
		edge [ source 1 target 4 label "-" ]
		edge [ source 1 target 6 label "-" ]
		edge [ source 2 target 3 label "=" ]
	]
	constrainAdj [ id 1 op "=" count 0 
		nodeLabels [ label "O" label "N" label "S" ]
		edgeLabels [ label "-" ]
	]
]
""")]