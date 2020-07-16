include("common.py")

def michaelAddition34Gen():
	r = RuleGen("Michael Addition ")
	r.left.extend([
		'# RC=C-',
		'edge [ source 1 target 3 label "=" ]',
	])
	r.context.extend([
		'# RC=C-',
		'node [ id 0 label "H" ]',
		'node [ id 1 label "C" ]',
		'node [ id 3 label "C" ]',
		'node [ id 4 label "H" ]',
		'node [ id 5 label "C" ]',
		'edge [ source 0 target 1 label "-" ]',
		'edge [ source 3 target 4 label "-" ]',
		'edge [ source 3 target 5 label "-" ]',
	])
	r.right.extend([
		'# RC=C-',
		'edge [ source 1 target 3 label "-" ]',
	])

	rBase = r
	for i in range(1, 5):
		rNum = rBase.clone()
		rNum.name += " %d" % i
		for r in attach_H_C(rNum, 1, 2):
			if i == 1 or i == 3:
				r.context.extend([
					'node [ id 6 label "O" ]',
					'edge [ source 5 target 6 label "=" ]'
				])
			elif i == 2 or i == 4:
				r.context.extend([
					'node [ id 6 label "N" ]',
					'edge [ source 5 target 6 label "#" ]'
				])
			else: assert False
			if i < 3:
				r.left.extend([
					'# Nu',
					'edge [ source 100 target 200 label "-" ]',
				])
				r.context.extend([
					'# Nu',
					'node [ id 100 label "H" ]',
				])
				r.right.extend([
					'# Nu',
					'edge [ source 1 target 200 label "-" ]',
					'edge [ source 3 target 100 label "-" ]',
				])
				for r in attach_Nu(r, 200):
					yield r.loadRule()
			else:
				r.left.extend([
					'# -EWG',
					'edge [ source 102 target 103 label "-" ]',
				])
				r.context.extend([
					'# -EWG',
					'node [ id 100 label "O" ]',
					'node [ id 101 label "C" ]',
					'node [ id 102 label "C" ]',
					'node [ id 103 label "H" ]',
					'edge [ source 100 target 101 label "=" ]',
					'edge [ source 101 target 102 label "-" ]',
				])
				r.right.extend([
					'edge [ source 1 target 102 label "-" ]',
					'edge [ source 3 target 103 label "-" ]',
				])
				for r in attach_H_C(r, 102, 104):
					for r in attach_EWG(r, 102, 200):
						yield r.loadRule()
michaelAddition34 = [a for a in michaelAddition34Gen()]
