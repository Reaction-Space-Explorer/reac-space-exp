# I put the rules for both deamination and addition (backwards reaction) in this single file

deamination_vic = ruleGMLString("""rule [
	ruleID "Deamination of Vicinal Diamine"
	labelType "term"
	left [
		edge [ source 1 target 2 label "-" ]
		edge [ source 1 target 3 label "-" ]
		edge [ source 1 target 6 label "-" ]
		edge [ source 2 target 7 label "-" ]
	]
	context [
		node [ id 1 label "C" ]
		node [ id 2 label "C" ]
		node [ id 3 label "N" ]
		edge [ source 3 target 4 label "-" ]
		node [ id 4 label "*" ]
		edge [ source 3 target 5 label "-" ]
		node [ id 5 label "*" ]
		node [ id 6 label "H" ]
		node [ id 7 label "H" ]
		edge [ source 2 target 8 label "-" ]
		node [ id 8 label "N" ]
		edge [ source 8 target 9 label "-" ]
		node [ id 9 label "*" ]
		edge [ source 8 target 10 label "-" ]
		node [ id 10 label "*" ]
	]
	right [
		edge [ source 1 target 2 label "=" ]
		edge [ source 1 target 7 label "-" ]
		edge [ source 3 target 6 label "-" ]
	]
	# None of the R groups should be bonded to =O or =N (that would be an amide/amidine)
	constrainAdj [ id 4 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" ]
	]
	constrainAdj [ id 5 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" ]
	]
	constrainAdj [ id 9 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" ]
	]
	constrainAdj [ id 10 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" ]
	]
]""")


deamination_vic_inv = ruleGMLString("""rule [
	ruleID "Deamination of Vicinal Diamine, (inverse)"
	labelType "term"
	left [
		edge [ source 1 target 2 label "=" ]
		edge [ source 1 target 7 label "-" ]
		edge [ source 3 target 6 label "-" ]
	]
	context [
		node [ id 1 label "C" ]
		node [ id 2 label "C" ]
		node [ id 3 label "N" ]
		edge [ source 3 target 4 label "-" ]
		node [ id 4 label "*" ]
		edge [ source 3 target 5 label "-" ]
		node [ id 5 label "*" ]
		node [ id 6 label "H" ]
		node [ id 7 label "H" ]
		edge [ source 2 target 8 label "-" ]
		node [ id 8 label "N" ]
		edge [ source 8 target 9 label "-" ]
		node [ id 9 label "*" ]
		edge [ source 8 target 10 label "-" ]
		node [ id 10 label "*" ]
	]
	right [
		edge [ source 1 target 2 label "-" ]
		edge [ source 1 target 3 label "-" ]
		edge [ source 1 target 6 label "-" ]
		edge [ source 2 target 7 label "-" ]
	]
	# None of the R groups should be bonded to =O or =N (that would be an amide/amidine)
	constrainAdj [ id 4 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" ]
	]
	constrainAdj [ id 5 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" ]
	]
	constrainAdj [ id 9 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" ]
	]
	constrainAdj [ id 10 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" ]
	]
]""")


deamination_two = ruleGMLString("""rule [
	ruleID "Deamination 2"
	labelType "term"
	left [
		edge [ source 1 target 3 label "-" ]
		edge [ source 1 target 5 label "-" ]
		edge [ source 3 target 4 label "-" ]
		edge [ source 2 target 6 label "-" ]
	]
	context [
		node [ id 1 label "C" ]
		edge [ source 1 target 2 label "-" ]
		node [ id 2 label "C" ]
		node [ id 3 label "O" ]
		node [ id 4 label "H" ]
		node [ id 5 label "H" ]
		node [ id 6 label "N" ]
		edge [ source 6 target 7 label "-" ]
		node [ id 7 label "*" ]
		edge [ source 6 target 8 label "-" ]
		node [ id 8 label "*" ]
	]
	right [
		edge [ source 1 target 3 label "=" ]
		edge [ source 4 target 6 label "-" ]
		edge [ source 2 target 5 label "-" ]
	]
	constrainAdj [ id 7 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" ]
	]
	constrainAdj [ id 8 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" ]
	]
]""")

deamination_3 = ruleGMLString("""rule [
	ruleID "Deamination 3"
	labelType "term"
	left [
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 3 target 4 label "-" ]
	]
	context [
		node [ id 1 label "H" ]
		node [ id 2 label "C" ]
		node [ id 3 label "C" ]
		node [ id 4 label "N" ]
		edge [ source 4 target 5 label "-" ]
		node [ id 5 label "*" ]
		edge [ source 4 target 6 label "-" ]
		node [ id 6 label "*" ]
	]
	right [
		edge [ source 2 target 3 label "=" ]
		edge [ source 1 target 4 label "-" ]
	]
	# C2 and C3 should not be bonded with O or N
	constrainAdj [ id 2 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "-" label "=" ]
	]
	constrainAdj [ id 3 op "=" count 1 
		nodeLabels [ label "O" label "N"]
		edgeLabels [ label "-" label "=" ]
	]
	# R groups should not be part of =O or =N (no amides or amidines should participate)
	constrainAdj [ id 5 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" ]
	]
	constrainAdj [ id 6 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" ]
	]
]""")

deamination_3_inv = ruleGMLString("""rule [
	ruleID "Deamination 3 (inverse)"
	labelType "term"
	left [
		edge [ source 2 target 3 label "=" ]
		edge [ source 1 target 4 label "-" ]
	]
	context [
		node [ id 1 label "H" ]
		node [ id 2 label "C" ]
		node [ id 3 label "C" ]
		node [ id 4 label "N" ]
		edge [ source 4 target 5 label "-" ]
		node [ id 5 label "*" ]
		edge [ source 4 target 6 label "-" ]
		node [ id 6 label "*" ]
	]
	right [
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 3 target 4 label "-" ]
	]
	# C2 and C3 should not be bonded with O or N
	constrainAdj [ id 2 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "-" label "=" ]
	]
	constrainAdj [ id 3 op "=" count 0 
		nodeLabels [ label "O" label "N"]
		edgeLabels [ label "-" label "=" ]
	]
	# R groups should not be part of =O or =N (no amides or amidines should participate)
	constrainAdj [ id 5 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" ]
	]
	constrainAdj [ id 6 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" ]
	]
]""")