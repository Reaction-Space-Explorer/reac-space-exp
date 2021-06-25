##
termo_oxigen3 = ruleGMLString("""rule [
	ruleID "termo O attack by .OH "  # O + .OH + Md -> .O-O-H + Md
	left [

		node [ id 1 label "O" ]
		node [ id 2 label "O." ]

		
		 
	]   
	context [
		node [ id 3 label "H" ]
		node [ id 4 label "Md" ]       		
      		edge [ source 2 target 3 label "-" ]		    		
		
	]   
	right [
		node [ id 1 label "O." ]
		node [ id 2 label "O" ]
		edge [ source 1 target 2 label "-" ]	
	
	]
	constrainAdj [ id 1 op "=" count 0 ]
		
		  
]""")

