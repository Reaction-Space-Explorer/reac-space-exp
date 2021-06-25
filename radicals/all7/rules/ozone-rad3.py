## For !) 
ozone3 = ruleGMLString("""rule [
	ruleID "O3 attack by .H "
	left [
		node [ id 1 label "O+" ]
		node [ id 2 label "O-" ]
		edge [ source 1 target 2 label "-" ]
		node [ id 3 label "H." ]
	]   
	context [
      		node [ id 4 label "*" ]
      		edge [ source 1 target 4 label "*" ]
		
	]   
	right [
		node [ id 1 label "O" ]
		node [ id 2 label "O." ]
		node [ id 3 label "H" ]
		edge [ source 2 target 3 label "-" ]
	]
	
]""")

