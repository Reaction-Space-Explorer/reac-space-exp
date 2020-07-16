cyanateFromCyanogen = [ruleGMLString("""rule [
	ruleID "Cyanate from cyanogen"
	left [
		# NCCN
		edge [ source 2 target 3 label "-" ]
		edge [ source 3 target 4 label "#" ]
		# Water
		edge [ source 5 target 6 label "-" ]
		edge [ source 5 target 7 label "-" ]
	]
	context [
		# Atoms from NCCN
		node [ id 1 label "N" ]
		edge [ source 1 target 2 label "#" ]
		node [ id 2 label "C" ]
		node [ id 3 label "C" ]
		node [ id 4 label "N" ]
		# Atoms from Water
		node [ id 5 label "O" ]
		node [ id 6 label "H" ]
		node [ id 7 label "H" ]
	]
	right [
		# HCN
		edge [ source 2 target 6 label "-" ]
		# Cyanate
		edge [ source 3 target 4 label "=" ]
		edge [ source 3 target 5 label "=" ]
		edge [ source 4 target 7 label "-" ]
	]
]""")]
