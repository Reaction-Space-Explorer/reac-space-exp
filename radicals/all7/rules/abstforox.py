## For !) 
abstract = ruleGMLString("""rule [
	ruleID "O3 attack by O "
	left [
		node [ id 1 label "O+" ]
		node [ id 2 label "O-" ]
		edge [ source 1 target 2 label "-" ]

	]   
	context [
      		node [ id 4 label "*" ]
      		edge [ source 1 target 4 label "*" ]
		node [ id 3 label "O" ]
	]   
	right [
		node [ id 1 label "O" ]
		node [ id 2 label "O" ]
		edge [ source 2 target 3 label "=" ]
	]
	
	constrainAdj [ id 3 op "=" count 0 ]	  
]""")

