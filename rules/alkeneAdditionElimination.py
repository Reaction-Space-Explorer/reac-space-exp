include("common.py")

def alkeneAdditionEliminationGen():
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
alkeneAdditionElimination = [a for a in alkeneAdditionEliminationGen()]