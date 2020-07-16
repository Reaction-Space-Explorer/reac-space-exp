hcnToKetone = [ruleGMLString("""rule [
	ruleID "HCN to Ketone" 
	left [
		edge [ source 2 target 3 label "=" ]
		edge [ source 5 target 6 label "-" ]
	]   
	context [
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 4 label "-" ]
		edge [ source 6 target 7 label "#" ]
		node [ id 1 label "C" ]
		node [ id 2 label "C" ]
		node [ id 3 label "O" ]
		node [ id 4 label "C" ]
		node [ id 5 label "H" ]
		node [ id 6 label "C" ]		
		node [ id 7 label "N" ]
	]	
	right [
		edge [ source 2 target 3 label "-" ]
		edge [ source 3 target 5 label "-" ]
		edge [ source 2 target 6 label "-" ]
	]   
]""")]
