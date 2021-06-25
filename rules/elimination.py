# This converts  COCO to C(=O)C (it converts the enol to keto in the same step)
#  to avoid enols
elimination1 = [ruleGMLString("""rule [
	ruleID "Elimination + enol to keto"
	left [ 
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 1 target 4 label "-" ]
		edge [ source 5 target 6 label "-" ]
	]
	context [
		node [ id 1 label "C" ]
		node [ id 2 label "O" ]
		node [ id 3 label "H" ]
		node [ id 4 label "H" ]
		node [ id 5 label "C" ]
		node [ id 6 label "O" ]
		node [ id 7 label "H" ]
		edge [ source 1 target 5 label "-" ]
		edge [ source 6 target 7 label "-" ]
	]
	right [
		edge [ source 1 target 2 label "=" ]
		edge [ source 3 target 6 label "-" ] # make H2O
		edge [ source 4 target 5 label "-" ] # Add -H to complete valency 
	]
	# neither carbon should be a part of a carboxylic acid/thiocarboxylic acid or anything
	constrainAdj [ id 1 op "=" count 0 
		nodeLabels [ label "O" label "N" label "S"]
		edgeLabels [ label "=" ]
	]
	constrainAdj [ id 5 op "=" count 0
		nodeLabels [ label "O" label "N" label "S"]
		edgeLabels [ label "=" ]
	]
	# it should not form amides
	constrainAdj [ id 1 op "=" count 0 
		nodeLabels [ label "N" label "S" ]
		edgeLabels [ label "-" ]
	]
]""")]


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
		nodeLabels [ label "O" ]
		edgeLabels [ label "-" ]
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


dehydration_amine = ruleGMLString("""rule [
	ruleID "Dehydration of Amines"
	labelType "term"
	left [
		edge [ source 1 target 2 label "-" ]
		edge [ source 1 target 3 label "-" ]
		edge [ source 2 target 5 label "-" ]
	]
	context [
		node [ id 1 label "C" ]
		node [ id 2 label "C" ]
		node [ id 3 label "O" ]
		edge [ source 3 target 4 label "-" ]
		node [ id 4 label "H" ]
		node [ id 5 label "H" ]
		node [ id 6 label "N" ]
		edge [ source 2 target 6 label "-" ]
		node [ id 7 label "*" ]
		edge [ source 6 target 7 label "-" ]
		node [ id 8 label "*" ]
		edge [ source 6 target 8 label "-" ]
	]
	right [
		edge [ source 1 target 2 label "=" ]
		edge [ source 3 target 5 label "-" ]
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

hydration_amine = ruleGMLString("""rule [
	ruleID "Deydration of Amines (inverse)"
	labelType "term"
	left [
		edge [ source 1 target 2 label "=" ]
		edge [ source 3 target 5 label "-" ]
	]
	context [
		node [ id 1 label "C" ]
		node [ id 2 label "C" ]
		node [ id 3 label "O" ]
		edge [ source 3 target 4 label "-" ]
		node [ id 4 label "H" ]
		node [ id 5 label "H" ]
		node [ id 6 label "N" ]
		edge [ source 2 target 6 label "-" ]
		node [ id 7 label "*" ]
		edge [ source 6 target 7 label "-" ]
		node [ id 8 label "*" ]
		edge [ source 6 target 8 label "-" ]
	]
	right [
		edge [ source 1 target 2 label "-" ]
		edge [ source 1 target 3 label "-" ]
		edge [ source 2 target 5 label "-" ]
	]
	# avoid making aminols
	constrainAdj [ id 1 op "=" count 0 
		nodeLabels [ label "N" ]
		edgeLabels [ label "-" ]
	] 
	constrainAdj [ id 2 op "=" count 0
		nodeLabels [ label "O" ]
		edgeLabels [ label "-" ]
	]

	# not amide/amidine, etc.
	constrainAdj [ id 7 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" ]
	]
	constrainAdj [ id 8 op "=" count 0 
		nodeLabels [ label "O" label "N" ]
		edgeLabels [ label "=" ]
	]
]""")

#hydration_amine.print()