include("common.py")
#A realizar una mejor busqueda de terminaciones, ver tambien la auto reaccion entre radicales 
def terminationGen():
	def terminCC(r, root, offset):
		
		#For CC 
		rCopy = r.clone()
		rCopy.name += ",CC"
		
		rCopy.left.extend([
		'node [ id %d label "C." ]' % (root),
		'node [ id %d label "C." ]' % (offset),

		])
		rCopy.right.extend([
		'node [ id %d label "C" ]' % (root),
		'node [ id %d label "C" ]' % (offset),
		'edge [ source %d target %d label "-" ]' % (root, offset),
		])
			
		yield rCopy

		#For CN
		rCopy = r.clone()
		rCopy.name += ",CN"
		
		rCopy.left.extend([
		'node [ id %d label "C." ]' % (root),
		'node [ id %d label "N." ]' % (offset),

		])
		rCopy.right.extend([
		'node [ id %d label "C" ]' % (root),
		'node [ id %d label "N" ]' % (offset),
		'edge [ source %d target %d label "-" ]' % (root, offset),
		])
			
		yield rCopy		
		
		
		#For CS
		rCopy = r.clone()
		rCopy.name += ",CS"
		
		rCopy.left.extend([
		'node [ id %d label "C." ]' % (root),
		'node [ id %d label "S." ]' % (offset),

		])
		rCopy.right.extend([
		'node [ id %d label "C" ]' % (root),
		'node [ id %d label "S" ]' % (offset),
		'edge [ source %d target %d label "-" ]' % (root, offset),
		])
			
		yield rCopy		
		
		#For CO
		rCopy = r.clone()
		rCopy.name += ",CO"
		
		rCopy.left.extend([
		'node [ id %d label "C." ]' % (root),
		'node [ id %d label "O." ]' % (offset),

		])
		rCopy.right.extend([
		'node [ id %d label "C" ]' % (root),
		'node [ id %d label "O" ]' % (offset),
		'edge [ source %d target %d label "-" ]' % (root, offset),
		])
			
		yield rCopy	
		
		#For NN
		rCopy = r.clone()
		rCopy.name += ",NN"
		
		rCopy.left.extend([
		'node [ id %d label "N." ]' % (root),
		'node [ id %d label "N." ]' % (offset),

		])
		rCopy.right.extend([
		'node [ id %d label "N" ]' % (root),
		'node [ id %d label "N" ]' % (offset),
		'edge [ source %d target %d label "-" ]' % (root, offset),
		])
			
		yield rCopy					


		#For NS
		rCopy = r.clone()
		rCopy.name += ",NS"
		
		rCopy.left.extend([
		'node [ id %d label "N." ]' % (root),
		'node [ id %d label "S." ]' % (offset),

		])
		rCopy.right.extend([
		'node [ id %d label "N" ]' % (root),
		'node [ id %d label "S" ]' % (offset),
		'edge [ source %d target %d label "-" ]' % (root, offset),
		])
			
		yield rCopy					

		#For NO
		rCopy = r.clone()
		rCopy.name += ",NO"
		
		rCopy.left.extend([
		'node [ id %d label "N." ]' % (root),
		'node [ id %d label "O." ]' % (offset),

		])
		rCopy.right.extend([
		'node [ id %d label "N" ]' % (root),
		'node [ id %d label "O" ]' % (offset),
		'edge [ source %d target %d label "-" ]' % (root, offset),
		])
			
		yield rCopy	


		#For SS
		rCopy = r.clone()
		rCopy.name += ",SS"
		
		rCopy.left.extend([
		'node [ id %d label "S." ]' % (root),
		'node [ id %d label "S." ]' % (offset),

		])
		rCopy.right.extend([
		'node [ id %d label "S" ]' % (root),
		'node [ id %d label "S" ]' % (offset),
		'edge [ source %d target %d label "-" ]' % (root, offset),
		])
			
		yield rCopy	

		#For SO
		rCopy = r.clone()
		rCopy.name += ",SO"
		
		rCopy.left.extend([
		'node [ id %d label "S." ]' % (root),
		'node [ id %d label "O." ]' % (offset),

		])
		rCopy.right.extend([
		'node [ id %d label "S" ]' % (root),
		'node [ id %d label "O" ]' % (offset),
		'edge [ source %d target %d label "-" ]' % (root, offset),
		])
		
		yield rCopy	
				

		#For OO
		rCopy = r.clone()
		rCopy.name += ",OO"
		
		rCopy.left.extend([
		'node [ id %d label "O." ]' % (root),
		'node [ id %d label "O." ]' % (offset),

		])
		rCopy.right.extend([
		'node [ id %d label "O" ]' % (root),
		'node [ id %d label "O" ]' % (offset),
		'edge [ source %d target %d label "-" ]' % (root, offset),
		])
		
		yield rCopy	
								
	r = RuleGen("termination")
	for r in terminCC(r, 0, 1):
		yield r.loadRule()

terminationcc = [a for a in terminationGen()]
