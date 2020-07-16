include("common.py")

def alkyneAdditionGen():
	r = RuleGen("Alkyne Addition")
	r.left.extend([
		'edge [ source 0 target 1 label "#" ]',
		'edge [ source 2 target 3 label "-" ]',
	])
	r.context.extend([
		'node [ id 0 label "C" ]',
		'node [ id 1 label "C" ]',
		'node [ id 2 label "H" ]',
	])
	r.right.extend([
		'edge [ source 0 target 1 label "=" ]',
		'edge [ source 0 target 2 label "-" ]',
		'edge [ source 1 target 3 label "-" ]',
	])
	for r in attach_Nu(r, 3):
		yield r.loadRule()
alkyneAddition = [a for a in alkyneAdditionGen()]
