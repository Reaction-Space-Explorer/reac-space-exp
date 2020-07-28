include("common.py")

ester_hydrolysis = ruleGMLString("""rule [
	ruleID "Ester Hydrolysis"
	labelType "term"
	left [
		edge [ source 2 target 4 label "-" ]
		edge [ source 6 target 8 label "-" ]
	]
	context [
		# The ester
		node [ id 1 label "C" ]
		node [ id 2 label "C" ]
		node [ id 3 label "O" ]
		node [ id 4 label "O" ]
		node [ id 5 label "C" ]
		# H2O
		node [ id 6 label "O" ]
		node [ id 7 label "H" ]
		node [ id 8 label "H" ]
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "=" ]
		edge [ source 4 target 5 label "-" ]
		edge [ source 6 target 7 label "-" ]
	]
	right [
		edge [ source 2 target 6 label "-" ]
		edge [ source 4 target 8 label "-" ]
	]
	# solely to avoid this rule from being inverible for now
	constrainAdj [ id 2 op "=" count 1
		nodeLabels [ label "O"]
		edgeLabels [ label "-" ]
	]
]""")
	

# Temporarily turn off ester formation (happens less often in basic medium)
'''def esterFormationHydrolysisExchangeGen():
	r = RuleGen("Ester Formation")
	r.label = "term"
	r.left.extend([
		'# OC(A)OR',
		'edge [ source 0 target 3 label "-" ]',
		'# HNRR',
		'edge [ source 100 target 101 label "-" ]',
	])
	r.context.extend([
		'# OC(A)OR',
		'node [ id 0 label "C" ]',
		'node [ id 1 label "O" ]',
		'node [ id 2 label "*" ]',
		'node [ id 3 label "O" ]',
		'edge [ source 0 target 1 label "=" ]',
		'edge [ source 0 target 2 label "-" ]',
		'# HOR',
		'node [ id 100 label "O" ]',
		'node [ id 101 label "H" ]',
		'node [ id 102 label "*" ]',
		'edge [ source 100 target 102 label "-" ]',
	])
	r.right.extend([
		'# OCN',
		'edge [ source 0 target 100 label "-" ]',
		'# HOR',
		'edge [ source 3 target 101 label "-" ]',
	])
	r.constraints.extend([
		'# TODO: not carbonyl',
		'constrainAdj [',
		'\tid 102 count 0 op "="',
		'\tnodeLabels [ label "O" ]',
		'\tedgeLabels [ label "=" ]',
		']',
	])
	for r1 in attach_H_C(r, 3, 4):
		yield r1.loadRule()


esterFormationHydrolysisExchange = [a for a in esterFormationHydrolysisExchangeGen()]'''
