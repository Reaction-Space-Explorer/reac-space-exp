include("common.py")

amide_hydrolysis = ruleGMLString("""rule [
	ruleID "Amide Hydrolysis"
	labelType "term"

	left [
		edge [ source 2 target 4 label "-" ]
		edge [ source 7 target 8 label "-" ]
	]
	context [
		node [ id 1 label "*" ]
		edge [ source 1 target 2 label "-" ]
		node [ id 2 label "C" ]
		edge [ source 2 target 3 label "=" ]
		node [ id 3 label "O" ]
		node [ id 4 label "N" ]
		edge [ source 4 target 5 label "-" ]
		node [ id 5 label "*" ]
		edge [ source 4 target 6 label "-" ]
		node [ id 6 label "*" ]

		node [ id 7 label "H" ]
		node [ id 8 label "O" ]
		edge [ source 8 target 9 label "-" ]
		node [ id 9 label "H" ]
	]
	right [
		edge [ source 7 target 4 label "-" ]
		edge [ source 8 target 2 label "-" ]
	]
]""")

amide_hydrolysis.print()

# The following does only hydrolysis, not formation. Hydrolysis is not the preferred direction
# in basic medium
def amideFormationHydrolysisGen():
	r = RuleGen("Amide Formation Hydrolysis, C")
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
		'node [ id 4 label "C" ]'
		'edge [ source 0 target 1 label "=" ]',
		'edge [ source 0 target 2 label "-" ]',
		'edge [ source 3 target 4 label "-" ]',
		'# HNRR',
		'node [ id 100 label "N" ]',
		'node [ id 101 label "H" ]',
		'node [ id 102 label "*" ]',
		'node [ id 103 label "*" ]',
		'edge [ source 100 target 102 label "-" ]',
		'edge [ source 100 target 103 label "-" ]',
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
		'# TODO: not carbonyl',
		'constrainAdj [',
		'\tid 103 count 0 op "="',
		'\tnodeLabels [ label "O" ]',
		'\tedgeLabels [ label "=" ]',
		']',
	])
	yield r.loadRule()
	'''for r1 in attach_H_C(r, 3, 4):
		yield r1.loadRule()'''

amideFormationHydrolysis= [a for a in amideFormationHydrolysisGen()]

p = GraphPrinter()
p.withColour = True
p.withIndex = True
for r in amideFormationHydrolysis:
	r.print(p)