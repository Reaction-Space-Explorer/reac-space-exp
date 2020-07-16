include("common.py")

def ketoneAldehydeToThioketoneThialGen():
	r = RuleGen("Ketone/Aldehyde to Thioketone/Thial")
	r.left.extend([
		'edge [ source 2 target 4 label "=" ]',
		'edge [ source 5 target 6 label "-" ]',
		'edge [ source 6 target 7 label "-" ]'
	])
	r.context.extend([
		'node [ id 2 label "C" ]',
		'node [ id 4 label "O" ]',
		'node [ id 5 label "H" ]',
		'node [ id 6 label "S" ]',
		'node [ id 7 label "H" ]'
	])
	r.right.extend([
		'edge [ source 4 target 5 label "-" ]',
		'edge [ source 4 target 7 label "-" ]',
		'edge [ source 2 target 6 label "=" ]'
	])
	for r1 in attach_2_H_C(r, 2, 0):
		yield r1.loadRule()
etoneAldehydeToThioketoneThial = [a for a in ketoneAldehydeToThioketoneThialGen()]
