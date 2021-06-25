include("common.py")
#Reacion de OH con moleculas de O=M
def hydroxylOMGen():
	def ohom(r, root, offset):
		
		#For OH. + O=O -> OH-O.=O 
		rCopy = r.clone()
		rCopy.name += ", OH. + O=O"
		
		rCopy.left.extend([
		'node [ id %d label "O." ]' % (offset),
		'node [ id %d label "O" ]' % (offset+1),		
		
		])		

		rCopy.context.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "O" ]' % (offset+2),
		'edge [ source %d target %d label "-" ]' % (root, offset),		
		'edge [ source %d target %d label "=" ]' % (offset+1, offset+2),
		])
		
		rCopy.right.extend([
		'node [ id %d label "O" ]' % (offset),
		'node [ id %d label "O." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		


		])
			
		yield rCopy

		#For OH. + C=O -> OH-C.=O 
		rCopy = r.clone()
		rCopy.name += ", OH. + C=O"
		
		rCopy.left.extend([
		'node [ id %d label "O." ]' % (offset),
		'node [ id %d label "C" ]' % (offset+1),		
		
		])		

		rCopy.context.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "O" ]' % (offset+2),
		'edge [ source %d target %d label "-" ]' % (root, offset),		
		'edge [ source %d target %d label "=" ]' % (offset+1, offset+2),
		])
		
		rCopy.right.extend([
		'node [ id %d label "O" ]' % (offset),
		'node [ id %d label "C." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		


		])
			
		yield rCopy



		#For OH. + N=O -> OH-N.=O 
		rCopy = r.clone()
		rCopy.name += ", OH. + N=O"
		
		rCopy.left.extend([
		'node [ id %d label "O." ]' % (offset),
		'node [ id %d label "N" ]' % (offset+1),		
		
		])		

		rCopy.context.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "O" ]' % (offset+2),
		'edge [ source %d target %d label "-" ]' % (root, offset),		
		'edge [ source %d target %d label "=" ]' % (offset+1, offset+2),
		])
		
		rCopy.right.extend([
		'node [ id %d label "O" ]' % (offset),
		'node [ id %d label "N." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		


		])
			
		yield rCopy

		#For OH. + S=O -> OH-S.=O 
		rCopy = r.clone()
		rCopy.name += ", OH. + S=O"
		
		rCopy.left.extend([
		'node [ id %d label "O." ]' % (offset),
		'node [ id %d label "S" ]' % (offset+1),		
		
		])		

		rCopy.context.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "O" ]' % (offset+2),
		'edge [ source %d target %d label "-" ]' % (root, offset),		
		'edge [ source %d target %d label "=" ]' % (offset+1, offset+2),
		])
		
		rCopy.right.extend([
		'node [ id %d label "O" ]' % (offset),
		'node [ id %d label "S." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		


		])
			
		yield rCopy

	
	r = RuleGen("OHRadOM")
	for r in ohom(r, 0, 1):
		yield r.loadRule()

ohradom = [a for a in hydroxylOMGen()]
