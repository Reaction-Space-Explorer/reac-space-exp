# Something to reconsider -- apaprently if we allow -NH2 to be attached to the carbons,
# we might end up getting too many products (when we detail with N containing comps)

elimination2 = [ruleGMLString("""rule [
	ruleID "Elimination2"
	left [
		edge [ source 1 target 2 label "-" ]
		edge [ source 3 target 4 label "-" ]
		edge [ source 1 target 3 label "-" ]
	]
	context [
		node [ id 1 label "C" ]
		node [ id 2 label "H" ]
		node [ id 3 label "C" ]
		node [ id 4 label "O" ]
		node [ id 5 label "H" ]
		edge [ source 4 target 5 label "-" ]
	]
	right [
		edge [ source 1 target 3 label "=" ]
		edge [ source 2 target 4 label "-" ]
	]
	constrainAdj [ id 1 op "=" count 0
		# This one should not be double bonded in general
		edgeLabels [ label "=" ]
	]
	# Avoid carbon bonded to N or S
	constrainAdj [ id 3 op "=" count 0
		nodeLabels [ label "N" label "S" ]
		edgeLabels [ label "-" label "=" ]
	]
	# The carbon shouldn't be a part of a carboxylic acid
	constrainAdj [ id 3 op "=" count 0
		nodeLabels [ label "O" ]
		edgeLabels [ label "=" ]
	]
]
""")]


'''betaElimination2 = [ruleGMLString("""
rule [
    ruleID "Beta Elimination2"
    left  [  
		edge [ source 1 target 4 label "-" ]
		edge [ source 1 target 5 label "-" ]
		edge [ source 8 target 5 label "-" ]
		edge [ source 5 target 6 label "-" ]
		edge [ source 9 target 10 label "-" ]
    ]
    context [
        node [ id 1 label "C" ]
		node [ id 2 label "O" ]
		node [ id 3 label "H" ]
		node [ id 4 label "H" ]
		node [ id 5 label "C" ]
		node [ id 6 label "O" ]
		node [ id 7 label "H" ]
        node [ id 8 label "H" ]
		node [ id 9 label "C" ]
		node [ id 10 label "O" ]
        node [ id 11 label "H" ]
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 6 target 7 label "-" ]
		edge [ source 9 target 5 label "-" ]
		edge [ source 10 target 11 label "-" ]
    ]
    right [
		edge [ source 6 target 4 label "-" ]
		edge [ source 1 target 5 label "=" ]
		edge [ source 8 target 5 label "-" ]
		edge [ source 9 target 10 label "-" ]
	]

]""")]'''