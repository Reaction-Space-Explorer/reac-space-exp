##
foto_2 = ruleGMLString("""rule [
	ruleID "foto O2 "  # H-O-H + Hf -> H-O. + H.
	left [

		node [ id 2 label "O" ]
		node [ id 3 label "H" ]	
		edge [source 2 target 3 label "-"]
		node [ id 4 label "Hf" ]  		 
	]   
	context [
		node [ id 1 label "H" ]       			
		edge [source 1 target 2 label "-"]
	]   
	right [
		node [ id 2 label "O." ]
		node [ id 3 label "H." ]	
	
	]
		
		  
]""")
