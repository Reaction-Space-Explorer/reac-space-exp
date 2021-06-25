##
foto_2 = ruleGMLString("""rule [
	ruleID "foto O3 "  # [O-]-[O+]=O + Hf -> O + O=O 
	left [

		node [ id 1 label "O-" ]
		node [ id 2 label "O+" ]	
		edge [source 1 target 2 label "-"]
		node [ id 4 label "Hf" ] 
	]   
	context [
		node [ id 3 label "O" ]		       			
		edge [source 2 target 3 label "="]
	]   
	right [
		node [ id 1 label "O" ]
		node [ id 2 label "O" ]		
	
	]
		
		  
]""")
