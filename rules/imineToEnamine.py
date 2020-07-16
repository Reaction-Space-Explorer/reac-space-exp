imineToEnamine = [ruleGMLString("""rule [
	ruleID "Imine to Enamine" 
    labelType "term"
	left [
		edge [ source 2 target 3 label "=" ]
		edge [ source 3 target 4 label "-" ]
		edge [ source 4 target 5 label "-" ]
	]   
	context [
		edge [ source 1 target 2 label "-" ]
		node [ id 1 label "*" ]
		node [ id 2 label "N" ]
		node [ id 3 label "C" ]
		node [ id 4 label "C" ]
		node [ id 5 label "H" ]
	]	
	right [
		edge [ source 2 target 5 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 3 target 4 label "=" ]
	]   
]""")]
