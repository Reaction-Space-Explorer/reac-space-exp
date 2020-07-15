##
termo_h2 = ruleGMLString("""rule [
	ruleID "termo .OH attack by .H  "  # H-O. + .H +Md -> H-O-H +Md
	left [

		node [ id 2 label "O." ]
		node [ id 3 label "H." ]	
		
		 
	]   
	context [
		node [ id 1 label "H" ]
		node [ id 4 label "Md" ]       		
	    	edge [source 1 target 2 label "-"]	
		
	]   
	right [
		node [ id 2 label "O" ]	
		node [ id 3 label "H" ]
		edge [source 2 target 3 label "-"]
	
	]
		
		  
]""")
