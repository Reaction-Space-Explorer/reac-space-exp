# Note: this is specifically for 5, 6 and 7 membered (including oxygen) rings; such a huge constraint
# has been put to reduce the number of products

def create_hemiacetal_formation(ring_size):
	rule = RuleGen(f"Hemiacetal Formation for {ring_size} membered rings")
	rule.left.extend([
		'edge [ source 2 target 3 label "-" ]',
		'edge [ source 4 target 5 label "=" ]'
		])
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
		print(f"Adding edges between {i} and {i+1}")
		rule.context.append(f'edge [ source {i} target {i+1} label "-" ]')
	print(f"adding edge between {ring_size+2} and 4")
	rule.context.append(f'edge [ source {ring_size+2} target 4 label "-" ]')
	rule.right.extend([
		'edge [ source 4 target 5 label "-" ]',
		'edge [ source 3 target 5 label "-" ]',
		'# close the ring by joining O and C',
		'edge [ source 4 target 2 label "-" ]'])

	return rule.loadRule()

hemiacetal_formation = [create_hemiacetal_formation(ring_size) for ring_size in (5, 6, 7)]