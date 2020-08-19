# Note: this is specifically for 5, 6 and 7 membered (including oxygen) rings; such a huge constraint
# has been put to reduce the number of products

def create_hemiacetal_formation(ring_size, is_inverse=False):
	"""
	A method for generalizing the hemiacetal formation rule to multiple ring sizes
	Keyword arguments:
	ring_size: number of atoms in the ring skeleton(e.g 5, 6 or 7)
	is_inverse: are we generating the inverse (backward direction) rule?

	Note: the second parameter is needed because rules with constraints can't be
	automatically inverted in the current version of MOD.

	return: the Rule object
	"""
	rule = RuleGen(f"Hemiacetal Formation for {ring_size} membered rings{', inverse' if is_inverse == True else ''}")
	if is_inverse == False:
		# The hemiacetal formation rule
		rule.left.extend([
			'edge [ source 2 target 3 label "-" ]',
			'edge [ source 4 target 5 label "=" ]'
			])
		rule.right.extend([
			'edge [ source 4 target 5 label "-" ]',
			'edge [ source 3 target 5 label "-" ]',
			'# close the ring by joining O and C',
			'edge [ source 4 target 2 label "-" ]'])
		# avoid applying if the C=O is part of an amide/carboxylic/thiocarboxylic acid
		rule.constraints.extend([
		'constrainAdj [ id 4 op "=" count 0',
		'\tnodeLabels [ label "O" label "N" label "S" ]',
		'\tedgeLabels [ label "-" ]'
		']'])
	else: # make inverse; left becomes right, right becomes left.
		rule.right.extend([
			'edge [ source 2 target 3 label "-" ]',
			'edge [ source 4 target 5 label "=" ]'
			])
		rule.left.extend([
			'edge [ source 4 target 5 label "-" ]',
			'edge [ source 3 target 5 label "-" ]',
			'edge [ source 4 target 2 label "-" ]'])
		# We need *some* constraints to make sure the rule doesn't invert itself
		# otherwise it would invert itself to the forward rule without constraints
		rule.constraints.extend([
		'constrainAdj [ id 4 op "<=" count 1',
		'\tnodeLabels [ label "O" label "N" label "S" ]',
		'\tedgeLabels [ label "-" ]'
		']'])
	
	# The general context
	rule.context.extend([
		'node [ id 1 label "C" ]',
		'node [ id 2 label "O" ]',
		'node [ id 3 label "H" ]',
		'node [ id 4 label "C" ]',
		'node [ id 5 label "O" ]',
		'edge [ source 1 target 2 label "-" ]'
		])
	for i in range(6, 6+ring_size-3): # 2 carbons already in, add n-3 more (one oxygen)
		rule.context.append(f'node [ id {i} label "C"] ')
	rule.context.append('edge [ source 1 target 6 label "-" ]')
	for i in range(6, 6+ring_size-4):
		rule.context.append(f'edge [ source {i} target {i+1} label "-" ]')
	rule.context.append(f'edge [ source {ring_size+2} target 4 label "-" ]')

	return rule.loadRule()

hemiacetal_formation = [create_hemiacetal_formation(ring_size) for ring_size in (5, 6, 7)]
hemiacetal_inverse = [create_hemiacetal_formation(ring_size, is_inverse=True)
						for ring_size in (5, 6, 7)]