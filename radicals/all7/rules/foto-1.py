##
foto_1 = ruleGMLString("""rule [
	ruleID "foto O2 "  # O=O + Hf -> O + O 
	left [

		node [ id 1 label "O" ]
		node [ id 2 label "O" ]	
		edge [source 1 target 2 label "="]
		node [ id 4 label "Hf" ]  
	]     
	right [
		node [ id 1 label "O" ]	
		node [ id 2 label "O" ]
	
	]
		
		  
]""")
