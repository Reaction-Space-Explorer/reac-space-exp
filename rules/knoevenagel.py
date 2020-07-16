include("common.py")

def knoevenagelGen():
	def attach_CN_CdOA(r, root, offset):
		r.context.extend([
			"# CN or C(=O)A",
			'edge [ source %d target %d label "-" ]' % (root, offset),
			'node [ id %d label "C" ]' % offset,
		])
		rCopy = r.clone()
		rCopy.name += ", CN"
		rCopy.context.extend([
			'node [ id %d label "N" ]' % (offset + 1),
			'edge [ source %d target %d label "#" ]' % (offset, offset + 1),
		])
		yield rCopy
		rCopy = r.clone()
		rCopy.name += ", C(=O)A"
		rCopy.label = "term"
		rCopy.context.extend([
			'node [ id %d label "O" ]' % (offset + 1),
			'edge [ source %d target %d label "=" ]' % (offset, offset + 1),
			'node [ id %d label "*" ]' % (offset + 2),
			'edge [ source %d target %d label "-" ]' % (offset, offset + 2),
		])
		yield rCopy


	r = RuleGen("Knoevenagel")
	r.left.extend([
		'# C-EWG',
		'edge [ source 0 target 1 label "-" ]',
		'# C(=O)',
		'edge [ source 10 target 11 label "=" ]',
	])
	r.context.extend([
		'# C-EWG',
		'node [ id 0 label "C" ]',
		'node [ id 1 label "H" ]',
		'# C=O',
		'node [ id 10 label "C" ]',
		'node [ id 11 label "O" ]',
	])
	r.right.extend([
		'edge [ source 0 target 10 label "-" ]',
		'edge [ source 1 target 11 label "-" ]',
		'edge [ source 10 target 11 label "-" ]',
	])
	# C-EWG
	for r in attach_H_C(r, 0, 100):
		for r in attach_CN_CdOA(r, 0, 200):
			for r in attach_EWG(r, 0, 300):
				# C=O
				r.name += " |"
				for r in attach_2_H_C(r, 10, 400):
					yield r.loadRule()

knoevenagel = [a for a in knoevenagelGen()]
