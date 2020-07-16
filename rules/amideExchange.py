amideExchange = [ruleGMLString("""rule [
	ruleID "Amide Exchange" 
    labelType "term"
	left [
		edge [ source 2 target 4 label "-" ]
		edge [ source 8 target 9 label "-" ]
	]   
	context [
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "=" ]
		edge [ source 4 target 5 label "-" ]
		edge [ source 4 target 6 label "-" ]
		edge [ source 7 target 8 label "-" ]
		edge [ source 8 target 10 label "-" ]
		node [ id 1 label "*" ]
		node [ id 2 label "C" ]
		node [ id 3 label "O" ]
		node [ id 4 label "N" ]
		node [ id 5 label "*" ]
		node [ id 6 label "*" ]
		node [ id 7 label "*" ]	
		node [ id 8 label "N" ]
		node [ id 9 label "H" ]
		node [ id 10 label "*" ]
	]
	right [
		edge [ source 2 target 8 label "-" ]
		edge [ source 4 target 9 label "-" ]
	]   
]""")]
