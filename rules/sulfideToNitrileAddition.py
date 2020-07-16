sulfideToNitrileAddition = [ruleGMLString("""rule [
	ruleID "Sulfide to Nitrile Addition" 
    labelType "term"
	left [
		edge [ source 2 target 3 label "#" ]
		edge [ source 4 target 5 label "-" ]
	] 
	context [
		edge [ source 1 target 2 label "-" ]
		edge [ source 5 target 6 label "-" ]
		node [ id 1 label "*" ]
		node [ id 2 label "C" ]
		node [ id 3 label "N" ]
		node [ id 4 label "H" ]
		node [ id 5 label "S" ]
		node [ id 6 label "*" ]
	]
	right [
		edge [ source 2 target 5 label "-" ]
		edge [ source 2 target 3 label "=" ]
		edge [ source 3 target 4 label "-" ]
	]   
]""")]
