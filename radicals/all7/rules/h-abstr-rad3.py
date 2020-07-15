##
habstr3 = ruleGMLString("""rule [
	ruleID ".OH attack by .O-X "
	left [
		node [ id 1 label "O." ]
		node [ id 3 label "O." ]
		edge [ source 1 target 2 label "-" ]
	]   
	context [
      		 node [ id 2 label "H" ]
      		 node [ id 4 label "*"]
		 edge [ source 3 target 4 label "-" ]
	]   
	right [
		node [ id 1 label "O" ]
		node [ id 3 label "O" ]
		edge [ source 2 target 3 label "-" ]

	]
	
		  
]""")

