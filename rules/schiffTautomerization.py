schiffTautomerization = [ruleGMLString("""rule [
	ruleID "Schiff Tautomerization" 
	left [
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 3 target 4 label "=" ]
	]   
	context [
		node [ id 1 label "H" ]
		node [ id 2 label "C" ]
		node [ id 3 label "N" ]
		node [ id 4 label "C" ]
	]	
	right [
		edge [ source 2 target 3 label "=" ]
		edge [ source 3 target 4 label "-" ]
		edge [ source 4 target 1 label "-" ]
	]   
]""")]
