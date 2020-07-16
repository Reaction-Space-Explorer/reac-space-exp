thioamidineToAmideHydrolysis = [ruleGMLString("""rule [
	ruleID "Thioamidine to Amide Hydrolysis" 
    labelType "term"
	left [
		edge [ source 2 target 3 label "-" ]
		edge [ source 2 target 5 label "=" ]
		edge [ source 7 target 8 label "-" ]
		edge [ source 8 target 9 label "-" ]
	]   
	context [
		
		edge [ source 1 target 2 label "-" ]
		edge [ source 3 target 4 label "-" ]
		edge [ source 5 target 6 label "-" ]
		node [ id 1 label "*" ]
		node [ id 2 label "C" ]
		node [ id 3 label "S" ]
		node [ id 4 label "*" ]
		node [ id 5 label "N" ]
		node [ id 6 label "*" ]
		node [ id 7 label "H" ]
		node [ id 8 label "O" ]
		node [ id 9 label "H" ]

	]	
	right [
		edge [ source 3 target 9 label "-" ]
		edge [ source 2 target 8 label "=" ]
		edge [ source 2 target 5 label "-" ]
		edge [ source 5 target 7 label "-" ]
	]   
]""")]
