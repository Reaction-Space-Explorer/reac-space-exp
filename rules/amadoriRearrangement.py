amadoriRearrangement = [ruleGMLString("""rule [
	ruleID "Amadori Rearrangement"
	left [
		# Big molecule
		edge [ source 2 target 5 label "-" ]
		# HNRR
		edge [ source 10 target 11 label "-" ]
	]
	context [
		# Big molecule
		node [ id 0 label "C" ]
		node [ id 1 label "C" ]
		node [ id 2 label "C" ]
		node [ id 3 label "C" ]
		node [ id 4 label "O" ]
		node [ id 5 label "O" ]
		node [ id 6 label "H" ]
		node [ id 7 label "H" ]
		edge [ source 0 target 1 label "-" ]
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 1 target 4 label "=" ]
		edge [ source 2 target 7 label "-" ]

		edge [ source 5 target 6 label "-" ]
		# HNRR
		node [ id 10 label "N" ]
		node [ id 11 label "H" ]
	]
	right [
		# Big molecule
		edge [ source 2 target 10 label "-" ]
		# Water
		edge [ source 5 target 11 label "-" ]
	]
]""")]
