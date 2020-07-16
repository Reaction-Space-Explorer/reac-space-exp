cannizarro2 = [ruleGMLString("""rule [
	ruleID "Cannizarro 2"
	labelType "term"
	left [
		edge [ source 8 target 7 label "-" ]
		edge [ source 2 target 3 label "=" ]
		edge [ source 6 target 10 label "-" ]
	]   
	context [
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 4 label "-" ]
		edge [ source 5 target 6 label "-" ]
		edge [ source 6 target 11 label "=" ]
		edge [ source 8 target 9 label "-" ]
		node [ id 1 label "*" ]
		node [ id 2 label "C" ]
		node [ id 3 label "O" ]
		node [ id 4 label "H" ]
		node [ id 5 label "*" ]
		node [ id 6 label "C" ]
		node [ id 7 label "H" ]
		node [ id 8 label "O" ]
		node [ id 9 label "H" ]
		node [ id 10 label "H" ]
		node [ id 11 label "O" ]
	]	
	right [
		edge [ source 6 target 8 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 3 target 10 label "-" ]
		edge [ source 2 target 7 label "-" ]
	]   
]""")]
