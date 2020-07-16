include("common.py")

def carbamylationGen():
	r = RuleGen("Carbamylation")
	r.left.extend([
		'# Isocyanate-H',
		'edge [ source 1 target 2 label "=" ]',
		'# The other',
		'edge [ source 105 target 150 label "-" ]',
	])
	r.context.extend([
		'# Isocyanate-H',
		'node [ id 0 label "O" ]',
		'node [ id 1 label "C" ]',
		'node [ id 2 label "N" ]',
		'node [ id 3 label "H" ]',
		'edge [ source 0 target 1 label "=" ]',
		'edge [ source 2 target 3 label "-" ]',
		'# The other',
		'node [ id 150 label "H" ]',
		'# NHR1 or SH',
	])
	r.right.extend([
		'edge [ source 1 target 2 label "-" ]',
		'edge [ source 2 target 150 label "-" ]',
		'edge [ source 1 target 105 label "-" ]',
	])
	rOld = r
	for atom in "SN":
		r = rOld.clone()
		r.name += ", %s" % atom
		r.context.extend([
			'node [ id 105 label "%s" ]' % atom,
		])
		if atom == "S":
			for end in "CH":
				r1 = r.clone()
				r1.name += end
				r1.context.extend([
					'# C or H',
					'node [ id 100 label "%s" ]' % end,
					'edge [ source 100 target 105 label "-" ]',
				])
				yield r1.loadRule()
			continue
		assert atom == "N"
		for r in attach_2_H_C(r, 105, 100, True):
			yield r.loadRule()

carbamylation = [a for a in carbamylationGen()]
