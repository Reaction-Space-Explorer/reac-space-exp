##
termo_h = ruleGMLString("""rule [
	ruleID "termo O2 attack by .H  "  # O=O + .H +Md -> .O-O-H +Md
	left [

		node [ id 1 label "O" ]
		node [ id 3 label "H." ]
      		edge [ source 1 target 2 label "=" ]	
		
		 
	]   
	context [
		node [ id 2 label "O" ]
		node [ id 4 label "Md" ]       		
	    		
		
	]   
	right [
		node [ id 1 label "O." ]
		node [ id 3 label "H" ]
      		edge [ source 1 target 2 label "-" ]
      		edge [ source 2 target 3 label "-" ]
	
	]
		
		  
]""")
