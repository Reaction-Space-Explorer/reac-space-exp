include("common.py")

def mannichGen():
	r = RuleGen("Mannich")
	r.left.extend([
		'# HC(=O)Ar/H',
		'edge [ source 0 target 1 label "=" ]',
		'# HN(R1)(R2)',
		'edge [ source 100 target 101 label "-" ]',
		'# O=C(R1)C(R2)(R3)H',
		'edge [ source 202 target 203 label "-" ]',
	])
	r.context.extend([
		'# HC(=O)Ar/H',
		'node [ id 0 label "C" ]',
		'node [ id 1 label "O" ]',
		'node [ id 2 label "H" ]',
		'edge [ source 0 target 2 label "-" ]',
		'edge [ source 0 target 3 label "-" ]',
		'# HN(R1)(R2)',
		'node [ id 100 label "N" ]',
		'node [ id 101 label "H" ]',
		'# O=C(R1)C(R2)(R3)H',
		'node [ id 200 label "O" ]',
		'node [ id 201 label "C" ]',
		'node [ id 202 label "C" ]',
		'node [ id 203 label "H" ]',
		'edge [ source 200 target 201 label "=" ]',
		'edge [ source 201 target 202 label "-" ]',
	])
	r.right.extend([
		'# Big molecule',
		'edge [ source 0 target 100 label "-" ]',
		'edge [ source 0 target 202 label "-" ]',
		'# water',
		'edge [ source 101 target 1 label "-" ]',
		'edge [ source 203 target 1 label "-" ]',
	])
	for a1 in ["H", "Ar"]:
		r1 = r.clone()
		r1.name += " %s" % a1
		if a1 == "Ar":
			r1.label = "term"
			r1.context.append('node [ id 3 label "*" ]')
			r1.constraints.extend([
				'constrainAdj [',
				'\tid 3 count 1 op ">="',
				'\tedgeLabels [ label ":" ]',
				']',
			])
		else:
			r1.context.append('node [ id 3 label "%s" ]' % a1)
		for r2 in attach_2_H_C(r1, 100, 102, True):
			for r3 in attach_H_C(r2, 201, 210):
				for r4 in attach_2_H_C(r3, 202, 220):
					yield r4.loadRule()

mannich = [a for a in mannichGen()]
