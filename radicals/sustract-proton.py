include("common.py")

def sustracProtonGen():
	def sustract_H_to_CorOorNorS(r, root, offset):
		
		#For C 
		
		rCopy = r.clone()
		rCopy.name += ",C"
		
		rCopy.left.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "C" ]' % (offset),
		'edge [ source %d target %d label "-" ]' % (root, offset),
		])
		rCopy.right.extend([
		'node [ id %d label "C." ]' % (offset),
		'node [ id %d label "H." ]' % (root),
		])
			
		yield rCopy
		
		#For N
		rCopy = r.clone()
		rCopy.name += ",N"
		
		rCopy.left.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "N" ]' % (offset),
		'edge [ source %d target %d label "-" ]' % (root, offset),
		])
		rCopy.right.extend([
		'node [ id %d label "N." ]' % (offset),
		'node [ id %d label "H." ]' % (root),
		])
		yield rCopy
		
		#For S
		rCopy = r.clone()
		rCopy.name += ",S"
		
		rCopy.left.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "S" ]' % (offset),
		'edge [ source %d target %d label "-" ]' % (root, offset),
		])
		rCopy.right.extend([
		'node [ id %d label "S." ]' % (offset),
		'node [ id %d label "H." ]' % (root),
		])		
		yield rCopy
		

		#For O
		rCopy = r.clone()
		rCopy.name += ",O"
		
		rCopy.left.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "O" ]' % (offset),
		'edge [ source %d target %d label "-" ]' % (root, offset),
		])
		rCopy.right.extend([
		'node [ id %d label "O." ]' % (offset),
		'node [ id %d label "H." ]' % (root),
		])		
		yield rCopy		
				
	r = RuleGen("sustract-Proton")
	for r in sustract_H_to_CorOorNorS(r, 0, 1):
		yield r.loadRule()

sustractproton = [a for a in sustracProtonGen()]
