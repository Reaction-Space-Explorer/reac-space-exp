include("common.py")

def gen_closure_nitrile(ring_size):
	r = RuleGen(f"Nitrile Ring Closure {ring_size}")
	r.label = "term"
	r.left.extend([
		'edge [ source 1 target 2 label "#" ]',
		'edge [ source 3 target 4 label "-" ]',
		'edge [ source 3 target 5 label "-" ]'
	])
	r.context.extend([
		'node [ id 1 label "N" ]',
		'node [ id 2 label "C" ]',
		'node [ id 3 label "N" ]',
		'node [ id 4 label "H" ]',
		'node [ id 5 label "H" ]',
		'node [ id 6 label "*" ]',
		'node [ id 7 label "*" ]',
		'node [ id 8 label "*" ]',
		'edge [ source 3 target 6 label "*" ]',
		'edge [ source 6 target 7 label "*" ]',
		'edge [ source 7 target 8 label "*" ]'
	])
	r.constraints.extend([
		'constrainAdj [ id 6 op "=" count 0',
		'\tnodeLabels [ label "O" ]',
		'\tedgeLabels [ label "=" ]'
		']'
	])
	if ring_size == 6: # could generalize to larger values but don't need at this point
		r.context.extend([
			'node [ id 9 label "*" ]',
			'edge [ source 8 target 9 label "*" ]',
			'edge [ source 2 target 9 label "*" ]'
		])
	elif ring_size == 5:
		r.context.extend(['edge [ source 8 target 2 label "*" ]'])
	r.right.extend([
		'edge [ source 3 target 2 label "=" ]',
		'edge [ source 1 target 2 label "-" ]',
		'edge [ source 1 target 4 label "-" ]',
		'edge [ source 1 target 5 label "-" ]'
	])
	return r.loadRule()

ring_closure_nitrile = [gen_closure_nitrile(size) for size in (5,6)]


def create_ring_closure(ring_size, is_inverse=False):
	"""
	A method for generating the Rule objects for ring closure
	Keyword arguments:
	ring_size: number of atoms in the ring skeleton(e.g 5, 6 or 7)
	is_inverse: are we generating the inverse (backward direction) rule?

	Note: the second parameter is needed because rules with constraints can't be
	automatically inverted in the current version of MOD.

	return: the Rule object
	"""
	# a list of rules, each involving a different element (see below)
	rules = []

	# so at the ends of the chain, will be electronegative atoms attached
	# either as a part of a C=O group (as in amide/etc. groups) or more directly (as amines)
	for atom1 in "ONS":
		for atom2 in "ONS":
			rule = RuleGen(f"Ring Closure {ring_size} membered{', inverse' if is_inverse == True else ''} {atom1}, {atom2}")
			rule.label = "term"
			rule.context.extend([
				'node [ id 1 label "C" ]',
				f'node [ id 2 label "{atom1}" ]',
				'node [ id 3 label "*" ]',
				'edge [ source 2 target 3 label "-" ]',
				# stuff on the other end of the ring-forming chain
				f'node [ id 4 label "{atom2}" ]',
				'node [ id 5 label "H" ]',
				'node [ id 6 label "O" ]',
				'edge [ source 1 target 6 label "=" ]'
			])
			# add other atoms of the chain
			for i in range(7, 7+ring_size-2): # 2 carbons already in, add n-3 more (one oxygen)
				rule.context.append(f'node [ id {i} label "C"] ')
			# connect these carbons in the chain together
			rule.context.append('edge [ source 1 target 7 label "-" ]')
			for i in range(7, 7+ring_size-3):
				rule.context.append(f'edge [ source {i} target {i+1} label "-" ]')
			rule.context.append(f'edge [ source {7+ring_size-3} target 4 label "-" ]')

			# Now do the left [] and right []
			# ring closure case (forward direction)
			if is_inverse == False:
				rule.left.extend([
					'edge [ source 1 target 2 label "-" ]',
					'edge [ source 4 target 5 label "-" ]'
				])
				rule.right.extend([
					'edge [ source 2 target 5 label "-" ]',
					'edge [ source 1 target 4 label "-" ]'
				])
			else: # right becomes left, left becomes right; add some contraints
				rule.right.extend([
					'edge [ source 1 target 2 label "-" ]',
					'edge [ source 4 target 5 label "-" ]'
				])
				rule.left.extend([
					'edge [ source 2 target 5 label "-" ]',
					'edge [ source 1 target 4 label "-" ]'
				])
				rule.constraints.extend([
					f'constrainAdj [ id {ring_size+4} op "<=" count 1',
					f'nodeLabels [ label "{atom2}"]',
					'edgeLabels [ label "-" ]'
					']'
				])
			rules.append(rule)

	return [r.loadRule() for r in rules]

ring_closure = []
for ring_size in (5, 6, 7):
	for r in create_ring_closure(ring_size):
		ring_closure.append(r)

ring_closure_inv = []
for ring_size in (5, 6, 7):
	for r in create_ring_closure(ring_size, is_inverse=True):
		ring_closure_inv.append(r)

# other ring closure rules.

amine_nitrile_5 = ruleGMLString("""rule [
	ruleID "Amine Nitrile Ring Closure 5"
	labelType "term"
	left [
		edge [ source 2 target 3 label "-" ] # break N-H
		edge [ source 7 target 8 label "#" ] # break C#N
	]
	context [
		node [ id 1 label "*" ] # C/H only
		node [ id 2 label "N" ]
		edge [ source 1 target 2 label "-" ]
		node [ id 3 label "H" ]
		node [ id 4 label "C" ]
		edge [ source 2 target 4 label "-" ]
		node [ id 5 label "*" ] # Anything
		edge [ source 4 target 5 label "*" ] # any type of bond
		node [ id 6 label "*" ]
		edge [ source 5 target 6 label "*" ]
		node [ id 7 label "C" ]
		edge [ source 6 target 7 label "-" ]
		node [ id 8 label "N" ]
	]
	right [
		edge [ source 7 target 8 label "=" ] # =N
		edge [ source 3 target 8 label "-" ] # =N-H
		edge [ source 2 target 7 label "-" ] # C-N (closes the ring)
	]
	# what the following constraint does is, make sure the C/H isn't attached to any N or O
	# apart from the one N it is already bonded to (atom 2). TODO: See if it deviates from the behavior
	constrainAdj [ id 1 op "=" count 1
		edgeLabels [ label "-" label "=" ]
		nodeLabels [ label "N" label "O" ]
	]
]""")


amine_nitrile_6 = ruleGMLString("""rule [
	ruleID "Amine Nitrile Ring Closure 6"
	labelType "term"
	left [
		edge [ source 2 target 3 label "-" ] # break N-H
		edge [ source 8 target 9 label "#" ] # break C#N
	]
	context [
		node [ id 1 label "*" ] # C/H only
		node [ id 2 label "N" ]
		edge [ source 1 target 2 label "-" ]
		node [ id 3 label "H" ]
		node [ id 4 label "C" ]
		edge [ source 2 target 4 label "-" ]
		node [ id 5 label "*" ] # Anything
		edge [ source 4 target 5 label "*" ] # any type of bond
		node [ id 6 label "*" ]
		edge [ source 5 target 6 label "*" ]
		node [ id 7 label "*" ]
		edge [ source 6 target 7 label "*" ]
		node [ id 8 label "C" ]
		edge [ source 7  target 8 label "-" ]
		node [ id 9 label "N" ]
	]
	right [
		edge [ source 8 target 9 label "=" ] # =N
		edge [ source 3 target 9 label "-" ] # =N-H
		edge [ source 2 target 8 label "-" ] # C-N (closes the ring)
	]
	# what the following constraint does is, make sure the C/H isn't attached to any N or O
	# apart from the one N it is already bonded to (atom 2). TODO: See if it deviates from the behavior
	constrainAdj [ id 1 op "=" count 1
		edgeLabels [ label "-" label "=" ]
		nodeLabels [ label "N" label "O" ]
	]
]""")


amine_imine_5 = ruleGMLString("""rule [
	ruleID "Amine Imine Ring Closure 5"
	labelType "term"
	left [
		edge [ source 3 target 4 label "-" ]
		edge [ source 9 target 10 label "-" ]
	]
	context [
		node [ id 1 label "H" ]
		node [ id 2 label "H" ]
		node [ id 3 label "N" ]
		edge [ source 1 target 3 label "-" ]
		edge [ source 2 target 3 label "-" ]
		node [ id 4 label "C" ]
		node [ id 5 label "*" ]
		edge [ source 4 target 5 label "*" ]
		node [ id 6 label "*" ]
		edge [ source 5 target 6 label "*" ]
		node [ id 7 label "C" ]
		edge [ source 6 target 7 label "-" ]
		node [ id 8 label "*" ]
		edge [ source 7 target 8 label "*" ]
		node [ id 9 label "N" ]
		edge [ source 7 target 9 label "=" ]
		node [ id 10 label "H" ]
	]
	right [
		edge [ source 4 target 9 label "-" ]
		edge [ source 10 target 3 label "-" ]
	]
]""")


amine_imine_6 = ruleGMLString("""rule [
	ruleID "Amine Imine Ring Closure 6"
	labelType "term"
	left [
		edge [ source 3 target 4 label "-" ]
		edge [ source 10 target 11 label "-" ]
	]
	context [
		node [ id 1 label "H" ]
		node [ id 2 label "H" ]
		node [ id 3 label "N" ]
		edge [ source 1 target 3 label "-" ]
		edge [ source 2 target 3 label "-" ]
		node [ id 4 label "C" ]
		node [ id 5 label "*" ]
		edge [ source 4 target 5 label "*" ]
		node [ id 6 label "*" ]
		edge [ source 5 target 6 label "*" ]
		node [ id 7 label "*" ]
		edge [ source 6 target 7 label "*" ]

		node [ id 8 label "C" ]
		edge [ source 7 target 8 label "-" ]
		node [ id 9 label "*" ]
		edge [ source 8 target 9 label "*" ]
		node [ id 10 label "N" ]
		edge [ source 8 target 10 label "=" ]
		node [ id 11 label "H" ]
	]
	right [
		edge [ source 4 target 10 label "-" ]
		edge [ source 11 target 3 label "-" ]
	]
]""")


# This rule gets printed horribly, but I think it's correct.
amidine_nitrile_6 = ruleGMLString("""rule [
	ruleID "Amidine Nitrile Ring Closure 6"
	labelType "term"
	left [
		edge [ source 1 target 2 label "#" ] #break C#N
		# remove some hydrogens
		edge [ source 5 target 6 label "-" ]
		edge [ source 9 target 10 label "-" ]
		# shift double bond
		edge [ source 5 target 7 label "-" ] # break N-C (to be converted into N=C)
		edge [ source 7 target 9 label "=" ] # break N=C (to be conv into N-C)
	]
	context [
		node [ id 1 label "N" ]
		node [ id 2 label "C" ]
		node [ id 3 label "*" ]
		edge [ source 2 target 3 label "-" ]
		node [ id 4 label "*" ]
		edge [ source 3 target 4 label "*" ]
		node [ id 5 label "N" ]
		edge [ source 4 target 5 label "-" ]
		node [ id 6 label "H" ]
		node [ id 7 label "C" ]
		node [ id 8 label "*" ]
		edge [ source 7 target 8 label "-" ]
		node [ id 9 label "N" ]
		node [ id 10 label "H" ]
	]
	right [
		edge [ source 5 target 7 label "=" ] # C=N
		edge [ source 7 target 9 label "-" ] # N-C
		edge [ source 2 target 9 label "=" ] # C=N
		# NH2
		edge [ source 1 target 2 label "-" ] # C-N
		edge [ source 1 target 6 label "-" ]
		edge [ source 1 target 10 label "-" ]
	]
]""")

amide_nitrile_6 = ruleGMLString("""rule [
	ruleID "Amide Nitrile Ring Closure 6"
	labelType "term"
	left [
		edge [ source 1 target 2 label "#" ] # break C#N
		edge [ source 9 target 10 label "-" ] # break N-H
		edge [ source 9 target 11 label "-" ]
	]
	context [
		node [ id 1 label "N" ]
		node [ id 2 label "C" ]
		node [ id 3 label "*" ]
		edge [ source 2 target 3 label "-" ]
		node [ id 4 label "*" ]
		edge [ source 3 target 4 label "*" ]
		node [ id 5 label "N" ]
		edge [ source 4 target 5 label "-"]
		node [ id 6 label "*" ]
		edge [ source 5 target 6 label "-" ]
		node [ id 7 label "C" ]
		edge [ source 5 target 7 label "-" ]
		node [ id 8 label "O" ]
		edge [ source 7 target 8 label "=" ]
		node [ id 9 label "N" ]
		edge [ source 7 target 9 label "-" ]
		node [ id 10 label "H" ]
		node [ id 11 label "H" ]
	]
	right [
		edge [ source 1 target 2 label "-" ]
		edge [ source 9 target 2 label "=" ]
		edge [ source 1 target 10 label "-" ]
		edge [ source 1 target 11 label "-" ]
	]
]""")


# I wrote it in a way that it will work for ring sizes 5-8
# but am currently using only 5 and 6
def create_amine_carbonyl(ring_size):
	# Need two rules, one with an O as a substituent and another with H.
	rule_h = RuleGen(f"Amine Carbonyl Ring Closure, H, {ring_size}")
	rule_c = RuleGen(f"Amine Carbonyl Ring Closure, C, {ring_size}")
	rule_c.label = "term"
	rule_h.label = "term"

	# most context will be common to both of these.
	comm_left = [
		'edge [ source 1 target 3 label "-" ]', # remove N-H
		'edge [ source 2 target 3 label "-" ]', # same as above
		'edge [ source 10 target 12 label "=" ]'
	]
	comm_context = [
		'node [ id 1 label "H" ]',
		'node [ id 2 label "H" ]',
		'node [ id 3 label "N" ]',
		'node [ id 4 label "C" ]',
		'edge [ source 3 target 4 label "-" ]',
		'node [ id 5 label "*" ]',
		'edge [ source 4 target 5 label "*" ]',
		'node [ id 6 label "*" ]',
		'edge [ source 5 target 6 label "*" ]',
		'node [ id 10 label "C" ]',
		f'edge [ source 10 target {ring_size+1} label "-" ]', # join '*' and C, this would depend on ring_size
		# atom 11 is C/H, defined separately for the two contexts
		'edge [ source 10 target 11 label "-" ]'
		'node [ id 12 label "O" ]'
	]
	comm_right = [
		'edge [ source 3 target 10 label "=" ]', # N=C
		'edge [ source 2 target 12 label "-" ]', # -H2O
		'edge [ source 1 target 12 label "-" ]' 
	]
	rule_c.left.extend(comm_left)
	rule_h.left.extend(comm_left)
	rule_c.right.extend(comm_right)
	rule_h.right.extend(comm_right)
	if ring_size > 5:
		for i in range(7, ring_size+2): # upto ring_size+1
			comm_context.extend([
				f'node [ id {i} label "*" ]',
				f'edge [ source {i-1} target {i} label "*" ]'
			])
	rule_c.context.extend(comm_context)
	rule_c.context.extend([
		'node [ id 11 label "C" ]'
	])
	rule_h.context.extend(comm_context)
	rule_h.context.extend([
		'node [ id 11 label "H" ]'
	])
	return [rule_h.loadRule(), rule_c.loadRule()]

amine_carbonyl_closure = []
# now generate the rules using the above method and add it to the list
for size in (5,6):
	rules = create_amine_carbonyl(size)
	for r in rules:
		amine_carbonyl_closure.append(r)


'''def ringClosureGen():
	def attachRingNode(r, id):
		rCopy = r.clone()
		rCopy.name += "C"
		rCopy.context.append('node [ id %d label "C" ]' % id)
		yield rCopy
	def attachRingBond(r, src, tar):
		rCopy = r.clone()
		rCopy.name += "-"
		rCopy.context.append('edge [ source %d target %d label "-" ]' % (src, tar))
		yield rCopy
	def endChain(rOrig, lastId):
		rOrig.context.append('node [ id %d label "C" ]' % lastId)
		for end in "SNO":
			r = rOrig.clone()
			r.left.append('edge [ source 200 target 201 label "-" ]')
			r.context.extend([
				'edge [ source %d target 200 label "-" ]' % lastId,
				'node [ id 200 label "%s" ]' % end,
				'node [ id 201 label "H" ]',
			])
			r.right.append('edge [ source 100 target 200 label "-" ]')
			if end == "N":
				for r in attach_H_C(r, 200, 202):
					for r in attachCOR(r):
						yield r
			else:
				for r in attachCOR(r):
					yield r
	def attachCOR(rOrig):
		for end in "ONS":
			r = rOrig.clone()
			r.left.append('edge [ source 100 target 102 label "-" ]')
			r.context.append('node [ id 102 label "%s" ]' % end)
			r.right.append('edge [ source 201 target 102 label "-" ]')
			if end == "O" or end == "S":
				for r in attach_H_C(r, 102, 103):
					yield r
			elif end == "N":
				for r in attach_2_H_C(r, 102, 103):
					yield r
			else:
				assert False
		
	# C* end: 100, ...
	# chain: 0, .., 10
	# other end: 200, ...
	for i in [5, 6]:
		r = RuleGen("Ring Closure %d, " % i)
		r.context.extend([
			'# TODO: add constraints',
			'# CO end',
			'node [ id 100 label "C" ]',
			'edge [ source 100 target 101 label "=" ]',
			'node [ id 101 label "O" ]',
			'# chain',
			'edge [ source 100 target 0 label "-" ]',
		])
		for r in attachRingNode(r, 0):
			for r in attachRingBond(r, 0, 1):
				for r in attachRingNode(r, 1):
					for r in attachRingBond(r, 1, 2):
						if i == 5:
							for r in endChain(r, 2):
								yield r.loadRule()
						else: # 6
							for r in attachRingNode(r, 2):
								for r in attachRingBond(r, 2, 3):
									for r in endChain(r, 3):
										yield r.loadRule()
ringClosure = [a for a in ringClosureGen()]'''
