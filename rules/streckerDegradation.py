include("common.py")

# TODO: is this the right trace?

def streckerDegradationGen():
	def attach_2_H_C(r, root, offset, withCarbonylConstraint=False):
		r.context.extend([
			"# 2, H or C on %d" % root,
		])
		for s1, s2 in [("H", "H"), ("H", "C"), ("C", "C")]:
			rCopy = r.clone()
			rCopy.name += ", " + s1 + s2
			rCopy.context.extend([
				'node [ id %d label "%s" ]' % (offset, s1),
				'node [ id %d label "%s" ]' % (offset + 1, s2),
			])
			rCopy.left.extend([
				'edge [ source %d target %d label "-" ]' % (root, offset),
				'edge [ source %d target %d label "-" ]' % (root, offset + 1),
			])
			rCopy.right.extend([
				'edge [ source 0 target %d label "-" ]' % offset,
				'edge [ source 0 target %d label "-" ]' % (offset + 1),
			])
			if s1 == 'C':
				rCopy.constraints.extend([
					'constrainAdj [',
					'	id %d op "=" count 0' % offset,
					'	nodeLabels [ label "O" ]',
					'	edgeLabels [ label "=" ]',
					']',
				])
			if s2 == 'C':
				rCopy.constraints.extend([
					'constrainAdj [',
					'	id %d op "=" count 0' % (offset + 1),
					'	nodeLabels [ label "O" ]',
					'	edgeLabels [ label "=" ]',
					']',
				])
			yield rCopy
	r = RuleGen("Strecker Degradation")
	r.label = "term"
	r.left.extend([
		'# CN',
		'edge [ source 0 target 2 label "-" ]',
		'edge [ source 0 target 10 label "-" ]',
		'edge [ source 10 target 12 label "-" ]',
		'edge [ source 12 target 13 label "-" ]',
		'# C=O',
	])
	r.context.extend([
		'# CN',
		'node [ id 0 label "C" ]',
		'node [ id 1 label "H" ]',
		'node [ id 2 label "*" ]',
		'node [ id 3 label "N" ]',
		'node [ id 4 label "H" ]',
		'node [ id 5 label "*" ]',
		'node [ id 10 label "C" ]',
		'node [ id 11 label "O" ]',
		'node [ id 12 label "O" ]',
		'node [ id 13 label "H" ]',
		'edge [ source 0 target 1 label "-" ]',
		'edge [ source 0 target 3 label "-" ]',
		'edge [ source 3 target 4 label "-" ]',
		'edge [ source 3 target 5 label "-" ]',
		'edge [ source 10 target 11 label "=" ]',
		'# C=O',
		'node [ id 100 label "C" ]',
		'node [ id 101 label "O" ]',
		'edge [ source 100 target 101 label "=" ]',
	])
	r.right.extend([
		'# AC=O',
		'edge [ source 100 target 2 label "-" ]',
		'edge [ source 100 target 13 label "-" ]',
		'# CO2',
		'edge [ source 10 target 12 label "=" ]',
		'# NC',
	])
	for r1 in attach_2_H_C(r, 100, 102):
		yield r1.loadRule()

streckerDegradation = [a for a in streckerDegradationGen()]