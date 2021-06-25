##
term_1 = ruleGMLString("""rule [
	ruleID "termination .O2H attack by H. " #H. + .O-O-H
	left [
		node [ id 1 label "H." ]
		node [ id 2 label "O." ]


	]   
	context [
		node [ id 3 label "O" ]
		node [ id 4 label "H" ] 		
 		edge [ source 2 target 3 label "-" ]	
 		edge [ source 3 target 4 label "-" ]	
		
	]   
	right [
		node [ id 1 label "H" ]
		node [ id 2 label "O" ]
		edge [ source 1 target 2 label "-" ]	
	]

		  
]""")

