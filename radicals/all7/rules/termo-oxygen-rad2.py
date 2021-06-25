##
termo_oxigen2 = ruleGMLString("""rule [
	ruleID "termo O attack by O "  # O + O + Md -> O=O + Md
	left [

		node [ id 4 label "Md" ] 
	]   
	context [
		node [ id 1 label "O" ]
		node [ id 2 label "O" ]
      		   		    		
		
	]   
	right [
		
		node [ id 4 label "Md" ] 
		edge [ source 1 target 2 label "=" ]
	
	]
	constrainAdj [ id 1 op "=" count 0 ]
	constrainAdj [ id 2 op "=" count 0 ]
		
		  
]""")

