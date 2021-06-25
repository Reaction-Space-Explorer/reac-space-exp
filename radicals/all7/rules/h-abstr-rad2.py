##
habstr2 = ruleGMLString("""rule [
	ruleID ".O2H attack by .O-X "
	left [
		node [ id 2 label "H" ]
		edge [ source 1 target 2 label "-" ]
		node [ id 3 label "O." ]
		node [ id 4 label "O." ]
		edge [ source 1 target 4 label "-" ]
	]   
	context [
      		 
      		node [ id 1 label "O" ]
      		node [ id 5 label "*" ]     		
		edge [ source 3 target 5 label "-" ]
	]   
	right [
		
		node [ id 2 label "H" ]
		edge [ source 2 target 3 label "-" ]
		node [ id 3 label "O" ]
		node [ id 4 label "O" ]
		edge [ source 1 target 4 label "=" ]
	]
	
		  
]""")

