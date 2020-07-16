include("common.py")

def cyanamidationGen():
	r = RuleGen("Cyanamidation")
	r.left.extend([
		'# H2NCN',
		'edge [ source 0 target 1 label "#" ]',
		'# The other',
		'edge [ source 105 target 150 label "-" ]',
	])
	r.context.extend([
		'# H2NCN',
		'node [ id 0 label "N" ]',
		'node [ id 1 label "C" ]',
		'node [ id 2 label "N" ]',
		'node [ id 3 label "H" ]',
		'node [ id 4 label "H" ]',
		'edge [ source 2 target 3 label "-" ]',
		'edge [ source 2 target 4 label "-" ]',
		'edge [ source 1 target 2 label "-" ]',
		'# The other',
		'node [ id 150 label "H" ]',
	])
	r.right.extend([
		'edge [ source 0 target 1 label "=" ]',
		'edge [ source 0 target 150 label "-" ]',
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

cyanamidation = [a for a in cyanamidationGen()]
