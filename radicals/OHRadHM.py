include("common.py")
#Reacion de OH con moleculas de HM
def hydroxylHMGen():
	def ohhm(r, root, offset):
		
		#For OH. + H-C 
		rCopy = r.clone()
		rCopy.name += "H-O. + H-C"
		
		rCopy.left.extend([
		'node [ id %d label "O." ]' % (offset),
		'node [ id %d label "C" ]' % (offset+2),		
		'edge [ source %d target %d label "-" ]' % (offset+1, offset+2),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "H" ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),		
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "O" ]' % (offset),
		'node [ id %d label "C." ]' % (offset+2),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		


		])
			
		yield rCopy

		#For OH. + H-N 
		rCopy = r.clone()
		rCopy.name += "H-O. + H-N"
		
		rCopy.left.extend([
		'node [ id %d label "O." ]' % (offset),
		'node [ id %d label "N" ]' % (offset+2),		
		'edge [ source %d target %d label "-" ]' % (offset+1, offset+2),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "H" ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),		
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "O" ]' % (offset),
		'node [ id %d label "N." ]' % (offset+2),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		


		])
			
		yield rCopy


		#For OH. + H-S 
		rCopy = r.clone()
		rCopy.name += "H-O. + H-S"
		
		rCopy.left.extend([
		'node [ id %d label "O." ]' % (offset),
		'node [ id %d label "S" ]' % (offset+2),		
		'edge [ source %d target %d label "-" ]' % (offset+1, offset+2),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "H" ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),		
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "O" ]' % (offset),
		'node [ id %d label "S." ]' % (offset+2),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		


		])
			
		yield rCopy

		#For OH. + H-O 
		rCopy = r.clone()
		rCopy.name += "H-O. + H-O"
		
		rCopy.left.extend([
		'node [ id %d label "O." ]' % (offset),
		'node [ id %d label "O" ]' % (offset+2),		
		'edge [ source %d target %d label "-" ]' % (offset+1, offset+2),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "H" ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),		
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "O" ]' % (offset),
		'node [ id %d label "O." ]' % (offset+2),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		


		])
			
		yield rCopy


				
	r = RuleGen("OHRadHM")
	for r in ohhm(r, 0, 1):
		yield r.loadRule()

ohradhm = [a for a in hydroxylHMGen()]
