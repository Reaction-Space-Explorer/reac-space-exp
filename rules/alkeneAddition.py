# This used to be an alkene addition/elimination rule file
# I replaced the O/N rules it had with customized ones
include("common.py")

# irreversible addition of HCN to C=C
hcn_addition_alk = ruleGMLString("""rule [
	ruleID "Alkene Addition, CN"
	labelType "term"
	left [
		edge [ source 1 target 2 label "=" ]
		edge [ source 3 target 4 label "-" ]
	]
	context [
		node [ id 1 label "C" ]
		node [ id 2 label "C" ]

		node [ id 3 label "H" ] # HCN
		node [ id 4 label "C" ]
		node [ id 5 label "N" ]
		edge [ source 4 target 5 label "#" ] 
	]
	right [
		edge [ source 1 target 2 label "-" ]
		edge [ source 1 target 3 label "-" ]
		edge [ source 2 target 4 label "-" ]
	]
]""")

'''def alkeneAdditionEliminationGen():
	r = RuleGen("Alkene Addition Elimination")
	r.left.extend([
		'edge [ source 0 target 1 label "-" ]',
		'edge [ source 1 target 2 label "-" ]',
		'edge [ source 2 target 3 label "-" ]',
	])
	r.context.extend([
		'node [ id 0 label "H" ]',
		'node [ id 1 label "C" ]',
		'node [ id 2 label "C" ]',
	])
	r.right.extend([
		'edge [ source 1 target 2 label "=" ]',
		'edge [ source 0 target 3 label "-" ]',
	])
	for r in attach_Nu(r, 3, allow_o=False):
		yield r.loadRule()
alkeneAdditionElimination = [a for a in alkeneAdditionEliminationGen()]'''