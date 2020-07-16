include("common.py")

def streckerDegradationDicarbonylGen():
	r = RuleGen("Strecker Degradation Dicarbonyl")
	r.left.extend([
		'edge [ source 1 target 2 label "-" ]',
		'edge [ source 2 target 3 label "-" ]',
		'edge [ source 3 target 5 label "-" ]',
		'edge [ source 5 target 8 label "-" ]',
		'edge [ source 14 target 15 label "=" ]',
		'edge [ source 12 target 14 label "-" ]'
	])
	r.context.extend([
		'edge [ source 3 target 4 label "=" ]',
		'edge [ source 5 target 6 label "-" ]',
		'edge [ source 8 target 9 label "-" ]',
		'edge [ source 12 target 13 label "=" ]',
		'node [ id 1 label "H" ]',
		'node [ id 2 label "O" ]',
		'node [ id 3 label "C" ]',
		'node [ id 4 label "O" ]',
		'node [ id 5 label "C" ]',
		'node [ id 6 label "H" ]',
		'node [ id 8 label "N" ]',
		'node [ id 9 label "H" ]',
		'node [ id 12 label "C" ]',
		'node [ id 13 label "O" ]',
		'node [ id 14 label "C" ]',
		'node [ id 15 label "O" ]',
	])
	r.right.extend([
		'edge [ source 14 target 1 label "-" ]',
		'edge [ source 12 target 14 label "-" ]',
		'edge [ source 5 target 15 label "=" ]',
		'edge [ source 2 target 3 label "=" ]',
		'edge [ source 14 target 8 label "-" ]'
	])
	for r1 in attach_H_C(r, 5, 7):
		for r2 in attach_H_C(r1, 8, 10):
			for r3 in attach_H_C(r2, 12, 11):
				for r4 in attach_H_C(r3, 14, 16):
					yield r4.loadRule()
streckerDegradationDicarbonyl = [a for a in streckerDegradationDicarbonylGen()]

