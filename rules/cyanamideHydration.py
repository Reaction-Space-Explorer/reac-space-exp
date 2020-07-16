include("common.py")

cyanamideHydration = [ruleGMLString("""rule [
	ruleID "Cyanamide Hydration"
	left [
		edge [ source 0 target 1 label "#" ]
		# Water
		edge [ source 10 target 11 label "-" ]
		edge [ source 10 target 12 label "-" ]
	]
	context [
		node [ id 0 label "N" ]
		node [ id 1 label "C" ]
		edge [ source 1 target 2 label "-" ]
		node [ id 2 label "N" ]
		edge [ source 2 target 3 label "-" ]
		node [ id 3 label "H" ]
		edge [ source 2 target 4 label "-" ]
		node [ id 4 label "H" ]

		node [ id 10 label "O" ]
		node [ id 11 label "H" ]
		node [ id 12 label "H" ]
	]
	right [
		edge [ source 0 target 1 label "-" ]
		edge [ source 0 target 11 label "-" ]
		edge [ source 0 target 12 label "-" ]

		edge [ source 1 target 10 label "=" ]
	]
]""")]
