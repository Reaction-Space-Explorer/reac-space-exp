include("common.py")

damnScission = [
	ruleGMLString("""rule [
	ruleID "DAMN Scission"
	left [
		edge [ source 0 target 10 label "=" ]	
		edge [ source 10 target 13 label "-" ]
		edge [ source 13 target 14 label "-" ]
		edge [ source 13 target 15 label "-" ]
	]
	context [
		# DAMN
		node [ id 0 label "C" ]

		edge [ source 0 target 1 label "-" ]
		node [ id 1 label "C" ]
		edge [ source 1 target 2 label "#" ]
		node [ id 2 label "N" ]

		edge [ source 0 target 3 label "-" ]
		node [ id 3 label "N" ]
		edge [ source 3 target 4 label "-" ]
		node [ id 4 label "H" ]
		edge [ source 3 target 5 label "-" ]
		node [ id 5 label "H" ]


		node [ id 10 label "C" ]

		edge [ source 10 target 11 label "-" ]
		node [ id 11 label "C" ]
		edge [ source 11 target 12 label "#" ]
		node [ id 12 label "N" ]

		node [ id 13 label "N" ]
		node [ id 14 label "H" ]
		node [ id 15 label "H" ]
	]
	right [
		edge [ source 10 target 13 label "#" ]
		edge [ source 0 target 14 label "-" ]
		edge [ source 0 target 15 label "-" ]
	]
]"""),
]