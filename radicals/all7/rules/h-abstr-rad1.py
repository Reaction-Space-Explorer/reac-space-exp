##
habstr1 = ruleGMLString("""rule [
	ruleID ".O2H attack by H. "
	left [
		node [ id 2 label "H" ]
		edge [ source 1 target 2 label "-" ]
		node [ id 3 label "H." ]
		node [ id 4 label "O." ]
		edge [ source 1 target 4 label "-" ]
	]   
	context [
      		 
      		node [ id 1 label "O" ]    		
		
	]   
	right [
		
		node [ id 2 label "H" ]
		edge [ source 2 target 3 label "-" ]
		node [ id 3 label "H" ]
		node [ id 4 label "O" ]
		edge [ source 1 target 4 label "=" ]
	]
	
	#constrainAdj [ id 3 op "=" count 0 ] 
		  
]""")

