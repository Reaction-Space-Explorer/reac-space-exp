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

]
""")]
# The following is Rana's elimination rule, the one I've written is quite different
'''betaElimination1 = [ruleGMLString("""
rule [
    ruleID "Elimination1"
    left  [  
		edge [ source 3 target 7 label "-" ]
		edge [ source 3 target 4 label "-" ]
		edge [ source 9 target 8 label "-" ]
		edge [ source 7 target 9 label "-" ]
    ]
    context [
        node [ id 1 label "C" ]
		node [ id 2 label "O" ]
        edge [ source 1 target 2 label "=" ]
		node [ id 3 label "C" ]
		edge [ source 1 target 3 label "-" ]
		node [ id 4 label "H" ]
		node [ id 7 label "C" ]
        node [ id 8 label "H" ]
		node [ id 9 label "O" ]
    ]
    right [
		edge [ source 3 target 7 label "=" ]
		edge [ source 9 target 4 label "-" ]
        edge [ source 9 target 8 label "-" ]
	]

]""")]
'''