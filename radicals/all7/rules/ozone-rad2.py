## For !) 
ozone2 = ruleGMLString("""rule [
	ruleID "O3 attack by .O-X, X != O "
	left [
		node [ id 1 label "O+" ]
		node [ id 2 label "O-" ]
		edge [ source 1 target 2 label "-" ]
		node [ id 3 label "O." ]
	]   
	context [
      		node [ id 4 label "*" ]
      		edge [ source 1 target 4 label "*" ]
		
	]   
	right [
		node [ id 1 label "O" ]
		node [ id 2 label "O." ]
		node [ id 3 label "O" ]
		edge [ source 2 target 3 label "-" ]
	]
	
	constrainAdj [id 3 op "=" count 0
		      nodeLabels [ label "O" ]
			]
]""")
