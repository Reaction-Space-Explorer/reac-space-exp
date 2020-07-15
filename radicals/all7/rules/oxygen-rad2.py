##
oxugen2 = ruleGMLString("""rule [
	ruleID "O attack by .O-OH " 
	left [
		node [ id 2 label "O." ]
		node [ id 3 label "O" ]
		edge [ source 2 target 3 label "-" ]
	]   
	context [
      		node [ id 1 label "O" ]     		    		
		
	]   
	right [
		
		node [ id 2 label "O" ]
		edge [ source 1 target 2 label "=" ]
		node [ id 3 label "O." ]
		
	]
	
	constrainAdj [ id 1 op "=" count 0 ] 
		  
]""")

