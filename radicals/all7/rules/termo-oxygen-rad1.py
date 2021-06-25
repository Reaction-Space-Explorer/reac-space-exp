##
termo_oxigen1 = ruleGMLString("""rule [
	ruleID "termo O=O attack by O "  # O + O=O + M -> [O-]-[O+]=O + M
	left [
		node [ id 1 label "O" ]
		node [ id 2 label "O" ]		
		node [ id 4 label "Md" ] 
	]   
	context [
      		node [ id 3 label "O" ]
      		edge [ source 2 target 3 label "=" ]   		    		
		
	]   
	right [
		
		node [ id 1 label "O-" ]
		node [ id 2 label "O+" ]
		node [ id 4 label "Md" ] 
		edge [ source 1 target 2 label "-" ]
	]
	
	constrainAdj [ id 1 op "=" count 0 ]	  
]""")

