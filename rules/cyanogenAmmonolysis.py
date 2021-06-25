include("common.py")

cyanogenAmonolysis = [ruleGMLString("""rule [
	ruleID "Cyanogen Ammonolysis"
	left [
		edge [ source 1 target 2 label "-" ]

		# Ammonium
		edge [ source 10 target 13 label "-" ]
	]
	context [
		node [ id 0 label "N" ]
		edge [ source 0 target 1 label "#" ]
		node [ id 1 label "C" ]

		node [ id 2 label "C" ]
		edge [ source 2 target 3 label "#" ]
		node [ id 3 label "N" ]


		node [ id 10 label "N" ]
		edge [ source 10 target 11 label "-" ]
		node [ id 11 label "H" ]
		edge [ source 10 target 12 label "-" ]
		node [ id 12 label "H" ]

		node [ id 13 label "H" ]
	]
	right [
		edge [ source 2 target 10 label "-" ]

		# HCN 
		edge [ source 1 target 13 label "-" ]
	]
]""")
]