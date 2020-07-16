ketoneCatalyzedDecarboxylation = [ruleGMLString("""rule [
	ruleID "Ketone-Catalyzed Decarboxylation"
    labelType "term"
	left [
		edge [ source 0 target 3 label "-" ]
		edge [ source 3 target 5 label "-" ]
		edge [ source 5 target 6 label "-" ]
	]
	context [
		node [ id 0 label "C" ]
		node [ id 1 label "*" ]
		node [ id 2 label "H" ]
		node [ id 3 label "C" ]
		node [ id 4 label "O" ]
		node [ id 5 label "O" ]
		node [ id 6 label "H" ]

		node [ id 7 label "N" ]
		node [ id 8 label "H" ]

		edge [ source 0 target 1 label "-" ]
		edge [ source 0 target 2 label "-" ]
		edge [ source 3 target 4 label "=" ]
		edge [ source 0 target 7 label "-" ]
		edge [ source 7 target 8 label "-" ]
	]
	right [
		edge [ source 0 target 6 label "-" ]
		edge [ source 3 target 5 label "=" ]
	]
]""")]
