include("common.py")

def gen_michael_addition(inverse=False):
	rules = []
	# The reaction is supposed to go S1 + S2 <=> P
	# define substrates as 
	substrate_1 = [
		[
			'node [ id 1 label "C" ]',
			'node [ id 2 label "C" ]',
			'node [ id 3 label "C" ]',
			'node [ id 4 label "O" ]',
			'edge [ source 2 target 3 label "-" ]',
			'edge [ source 3 target 4 label "=" ]'
		],
		[
			'node [ id 1 label "C" ]',
			'node [ id 2 label "C" ]',
			'node [ id 3 label "C" ]',
			'node [ id 4 label "O" ]',
			'edge [ source 2 target 3 label "-" ]',
			'edge [ source 3 target 4 label "=" ]'
		],
		[
			'node [ id 1 label "C" ]',
			'node [ id 2 label "C" ]',
			'node [ id 3 label "C" ]',
			'node [ id 4 label "N" ]',
			'edge [ source 2 target 3 label "-" ]',
			'edge [ source 3 target 4 label "#" ]'
		],
		[
			'node [ id 1 label "C" ]',
			'node [ id 2 label "C" ]',
			'node [ id 3 label "C" ]',
			'node [ id 4 label "N" ]',
			'edge [ source 2 target 3 label "-" ]',
			'edge [ source 3 target 4 label "#" ]'
		]
	]

	substrate_2 = [
		[
			'node [ id 11 label "H" ]',
			'node [ id 12 label "N" ]',
			'node [ id 13 label "*" ]',
			'node [ id 14 label "*" ]',
			'edge [ source 12 target 13 label "-" ]',
			'edge [ source 12 target 14 label "-" ]'
		],
		[
			'node [ id 11 label "H" ]',
			'node [ id 12 label "S" ]',
			'node [ id 13 label "*" ]',
			'edge [ source 12 target 13 label "-" ]'
		],
		[
			'node [ id 11 label "H" ]',
			'node [ id 12 label "C" ]',
			'node [ id 13 label "*" ]',
			'node [ id 14 label "*" ]',
			'edge [ source 12 target 13 label "-" ]',
			'edge [ source 12 target 14 label "-" ]'
		]
	]
	for s1 in range(len(substrate_1)):
		for s2 in range(len(substrate_2)):
			if s1 % 2 == 1 and inverse == True: # breaking of C#C should be irreversible
				break
			rule = RuleGen(f"Michael Addition {s1},{s2}, {'(reverse)' if inverse==True else ''}")
			rule.context.extend(substrate_1[s1])
			rule.context.extend(substrate_2[s2])
			if inverse == False:
				rule.left.extend([
					'edge [ source 11 target 12 label "-" ]'
				])
				# well in the substrate_1 list, item 1 has C=C, 2 has C#C, 3 has C=C and 4 again has C#C
				# so i need to treat them differently
				if s1 % 2 == 0: # these have C=C 
					rule.left.append('edge [ source 1 target 2 label "=" ]') # break the C=C
					rule.right.append('edge [ source 1 target 2 label "-" ]') # form C-C
				else:
					rule.left.append('edge [ source 1 target 2 label "#" ]') # break C#C
					rule.right.append('edge [ source 1 target 2 label "=" ]')
				rule.right.extend([
					'edge [ source 12 target 1 label "-" ]', # addition
					'edge [ source 11 target 2 label "-" ]' # connect the hydrogen
 				])
			else:
				rule.right.extend([
					'edge [ source 11 target 12 label "-" ]'
				])
				# well in the substrate_1 list, item 1 has C=C, 2 has C#C, 3 has C=C and 4 again has C#C
				# so i need to treat them differently
				if s1 % 2 == 0: # these have C=C 
					rule.left.append('edge [ source 1 target 2 label "-" ]')
					rule.right.append('edge [ source 1 target 2 label "=" ]')
				else:
					rule.left.append('edge [ source 1 target 2 label "=" ]')
					rule.right.append('edge [ source 1 target 2 label "#" ]')
				rule.left.extend([
					'edge [ source 12 target 1 label "-" ]', 
					'edge [ source 11 target 2 label "-" ]'
 				])
			# Add constraints
			if s1 % 2 == 0: # if there's a C=C formation involved, don't invert if enols are a possibility
				rule.constraints.extend([
					'constrainAdj [',
					'\tid 1 op "=" count 0',
					'\tnodeLabels [ label "O" ]',
					'\tedgeLabels [ label "-" ]',
					']',
					'constrainAdj [',
					'\tid 2 op "=" count 0',
					'\tnodeLabels [ label "O" ]',
					'\tedgeLabels [ label "-" ]',
					']'
				])
			if s2 == 0:
				rule.constraints.extend([
					'constrainAdj [',
					'\tid 1 op "=" count 0', # avoid forming diamine or aminols
					'\tnodeLabels [ label "O" label "N" ]',
					'\tedgeLabels [ label "-" ]',
					']',
					'constrainAdj [',
					'\tid 13 op "=" count 0',
					'\tnodeLabels [ label "O" ]',
					'\tedgeLabels [ label "=" ]',
					']',
					'constrainAdj [',
					'\tid 14 op "=" count 0',
					'\tnodeLabels [ label "O" ]',
					'\tedgeLabels [ label "=" ]',
					']',
					'constrainAdj [',
					'\tid 12 op "=" count 0',
					'\tnodeLabels [ label "O" label "N" label "S" ]',
					'\tedgeLabels [ label "-" ]',
					']'
				])
			elif s2 == 1:
				rule.constraints.extend([
					'constrainAdj [',
					'\tid 13 op "=" count 0',
					'\tnodeLabels [ label "O" ]',
					'\tedgeLabels [ label "=" ]',
					']',
				])
			elif s2 == 2:
				rule.constraints.extend([
					'constrainAdj [',
					'\tid 13 op "=" count 0',
					'\tnodeLabels [ label "O" label "N" ]',
					'\tedgeLabels [ label "=" label "#" ]',
					']',
				])
			rules.append(rule)
	return [r.loadRule() for r in rules]

michael_addition = gen_michael_addition()
michael_addition_inv = gen_michael_addition(inverse=True)

'''def michaelAddition34Gen():
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
michaelAddition34 = [a for a in michaelAddition34Gen()]'''