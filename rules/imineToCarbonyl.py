imineToCarbonyl = [ruleGMLString("""rule [
	ruleID "Imine to Carbonyl" 
    labelType "term"
	left [
		edge [ source 2 target 4 label "=" ]
		edge [ source 6 target 7 label "-" ]
		edge [ source 7 target 8 label "-" ]
	]   
	context [
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 4 target 5 label "-" ]
		node [ id 1 label "*" ]
		node [ id 2 label "C" ]
		node [ id 3 label "*" ]
		node [ id 4 label "N" ]
		node [ id 5 label "*" ]
		node [ id 6 label "H" ]
		node [ id 7 label "O" ]
		node [ id 8 label "H" ]
	]	
	right [
		edge [ source 2 target 7 label "=" ]
		edge [ source 4 target 6 label "-" ]
		edge [ source 4 target 8 label "-" ]
	]   
]""")]
