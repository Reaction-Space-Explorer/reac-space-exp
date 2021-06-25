##
oxugen1 = ruleGMLString("""rule [
	ruleID "O attack by .O-H " 
	left [
		node [ id 2 label "O." ]
		node [ id 3 label "H" ]
		edge [ source 2 target 3 label "-" ]
	]   
	context [
      		node [ id 1 label "O" ]     		    		
		
	]   
	right [
		
		node [ id 2 label "O" ]
		edge [ source 1 target 2 label "=" ]
		node [ id 3 label "H." ]
		
	]
	
	constrainAdj [ id 1 op "=" count 0 ] 
		  
]""")

