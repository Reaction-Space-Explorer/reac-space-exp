include("common.py")
#Reacion de un radical en un doble enlace
def doublebondattackGen():
	def doublemrad(r, root, offset):
		
		#For C. + C=C -> C-C-C.
		rCopy = r.clone()
		rCopy.name += "C. + C=C"
		
		rCopy.left.extend([
		'node [ id %d label "C." ]' % (root),
		'node [ id %d label "C" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "C" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "C" ]' % (root),
		'node [ id %d label "C." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy


		
		#For C. + C=O -> C-C-O.
		rCopy = r.clone()
		rCopy.name += "C. + C=O"
		
		rCopy.left.extend([
		'node [ id %d label "C." ]' % (root),
		'node [ id %d label "O" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "C" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "C" ]' % (root),
		'node [ id %d label "O." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy

		
		#For C. + C=N -> C-C-N.
		rCopy = r.clone()
		rCopy.name += "C. + C=N"
		
		rCopy.left.extend([
		'node [ id %d label "C." ]' % (root),
		'node [ id %d label "N" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "C" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "C" ]' % (root),
		'node [ id %d label "N." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy



		#For C. + C=S -> C-C-S.
		rCopy = r.clone()
		rCopy.name += "C. + C=S"
		
		rCopy.left.extend([
		'node [ id %d label "C." ]' % (root),
		'node [ id %d label "S" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "C" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "C" ]' % (root),
		'node [ id %d label "S." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy
		
##################


		
		#For C. + N=C -> C-N-C.
		rCopy = r.clone()
		rCopy.name += "C. + N=C"
		
		rCopy.left.extend([
		'node [ id %d label "C." ]' % (root),
		'node [ id %d label "C" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "N" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "C" ]' % (root),
		'node [ id %d label "C." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy


		
		#For C. + N=O -> C-N-O.
		rCopy = r.clone()
		rCopy.name += "C. + N=O"
		
		rCopy.left.extend([
		'node [ id %d label "C." ]' % (root),
		'node [ id %d label "O" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "N" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "C" ]' % (root),
		'node [ id %d label "O." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy


		
		#For C. + N=N -> C-N-N.
		rCopy = r.clone()
		rCopy.name += "C. + N=N"
		
		rCopy.left.extend([
		'node [ id %d label "C." ]' % (root),
		'node [ id %d label "N" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "N" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "C" ]' % (root),
		'node [ id %d label "N." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy


		
		#For C. + N=S -> C-N-S.
		rCopy = r.clone()
		rCopy.name += "C. + N=S"
		
		rCopy.left.extend([
		'node [ id %d label "C." ]' % (root),
		'node [ id %d label "S" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "N" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "C" ]' % (root),
		'node [ id %d label "S." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy
		
##################
		
		#For C. + S=C -> C-S-C.
		rCopy = r.clone()
		rCopy.name += "C. + S=C"
		
		rCopy.left.extend([
		'node [ id %d label "C." ]' % (root),
		'node [ id %d label "C" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "S" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "C" ]' % (root),
		'node [ id %d label "C." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy


		
		#For C. + S=O -> C-S-O.
		rCopy = r.clone()
		rCopy.name += "C. + S=O"
		
		rCopy.left.extend([
		'node [ id %d label "C." ]' % (root),
		'node [ id %d label "O" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "S" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "C" ]' % (root),
		'node [ id %d label "O." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy


		
		#For C. + S=N -> C-S-N.
		rCopy = r.clone()
		rCopy.name += "C. + S=N"
		
		rCopy.left.extend([
		'node [ id %d label "C." ]' % (root),
		'node [ id %d label "N" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "S" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "C" ]' % (root),
		'node [ id %d label "N." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy


		
		#For C. + S=S -> C-S-S.
		rCopy = r.clone()
		rCopy.name += "C. + S=S"
		
		rCopy.left.extend([
		'node [ id %d label "C." ]' % (root),
		'node [ id %d label "S" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "S" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "C" ]' % (root),
		'node [ id %d label "S." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy
		

#################

		#For C. + O=C -> C-O-C.
		rCopy = r.clone()
		rCopy.name += "C. + O=C"
		
		rCopy.left.extend([
		'node [ id %d label "C." ]' % (root),
		'node [ id %d label "C" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "O" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "C" ]' % (root),
		'node [ id %d label "C." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy


		
		#For C. + O=O -> C-O-O.
		rCopy = r.clone()
		rCopy.name += "C. + O=O"
		
		rCopy.left.extend([
		'node [ id %d label "C." ]' % (root),
		'node [ id %d label "O" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "O" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "C" ]' % (root),
		'node [ id %d label "O." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy


		
		#For C. + O=N -> C-O-N.
		rCopy = r.clone()
		rCopy.name += "C. + O=N"
		
		rCopy.left.extend([
		'node [ id %d label "C." ]' % (root),
		'node [ id %d label "N" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "O" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "C" ]' % (root),
		'node [ id %d label "N." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy


		
		#For C. + O=S -> C-O-S.
		rCopy = r.clone()
		rCopy.name += "C. + O=S"
		
		rCopy.left.extend([
		'node [ id %d label "C." ]' % (root),
		'node [ id %d label "S" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "O" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "C" ]' % (root),
		'node [ id %d label "S." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy
		
##############
##############

		
		#For H. + C=C -> C-C-C.
		rCopy = r.clone()
		rCopy.name += "H. + C=C"
		
		rCopy.left.extend([
		'node [ id %d label "H." ]' % (root),
		'node [ id %d label "C" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "C" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "C." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy


		
		#For H. + C=O -> C-C-O.
		rCopy = r.clone()
		rCopy.name += "H. + C=O"
		
		rCopy.left.extend([
		'node [ id %d label "H." ]' % (root),
		'node [ id %d label "O" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "C" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "O." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy

		
		#For H. + C=N -> H-C-N.
		rCopy = r.clone()
		rCopy.name += "H. + C=N"
		
		rCopy.left.extend([
		'node [ id %d label "H." ]' % (root),
		'node [ id %d label "N" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "C" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "N." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy



		#For H. + C=S -> H-C-S.
		rCopy = r.clone()
		rCopy.name += "H. + C=S"
		
		rCopy.left.extend([
		'node [ id %d label "H." ]' % (root),
		'node [ id %d label "S" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "C" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "S." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy
		
##################


		
		#For H. + N=C -> H-N-C.
		rCopy = r.clone()
		rCopy.name += "H. + N=C"
		
		rCopy.left.extend([
		'node [ id %d label "H." ]' % (root),
		'node [ id %d label "C" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "N" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "C." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy


		
		#For H. + N=O -> H-N-O.
		rCopy = r.clone()
		rCopy.name += "H. + N=O"
		
		rCopy.left.extend([
		'node [ id %d label "H." ]' % (root),
		'node [ id %d label "O" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "N" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "O." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy


		
		#For H. + N=N -> H-N-N.
		rCopy = r.clone()
		rCopy.name += "H. + N=N"
		
		rCopy.left.extend([
		'node [ id %d label "H." ]' % (root),
		'node [ id %d label "N" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "N" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "N." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy


		
		#For H. + N=S -> C-N-S.
		rCopy = r.clone()
		rCopy.name += "H. + N=S"
		
		rCopy.left.extend([
		'node [ id %d label "H." ]' % (root),
		'node [ id %d label "S" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "N" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "S." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy
		
##################
		
		#For H. + S=C -> H-S-C.
		rCopy = r.clone()
		rCopy.name += "H. + S=C"
		
		rCopy.left.extend([
		'node [ id %d label "H." ]' % (root),
		'node [ id %d label "C" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "S" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "C." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy


		
		#For H. + S=O -> H-S-O.
		rCopy = r.clone()
		rCopy.name += "H. + S=O"
		
		rCopy.left.extend([
		'node [ id %d label "H." ]' % (root),
		'node [ id %d label "O" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "S" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "O." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy


		
		#For H. + S=N -> H-S-N.
		rCopy = r.clone()
		rCopy.name += "H. + S=N"
		
		rCopy.left.extend([
		'node [ id %d label "H." ]' % (root),
		'node [ id %d label "N" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "S" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "N." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy


		
		#For H. + S=S -> H-S-S.
		rCopy = r.clone()
		rCopy.name += "H. + S=S"
		
		rCopy.left.extend([
		'node [ id %d label "H." ]' % (root),
		'node [ id %d label "S" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "S" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "S." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy
		

#################

		#For H. + O=C -> H-O-C.
		rCopy = r.clone()
		rCopy.name += "H. + O=C"
		
		rCopy.left.extend([
		'node [ id %d label "H." ]' % (root),
		'node [ id %d label "C" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "O" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "C." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy


		
		#For H. + O=O -> H-O-O.
		rCopy = r.clone()
		rCopy.name += "H. + O=O"
		
		rCopy.left.extend([
		'node [ id %d label "H." ]' % (root),
		'node [ id %d label "O" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "O" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "O." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy


		
		#For H. + O=N -> C-O-N.
		rCopy = r.clone()
		rCopy.name += "H. + O=N"
		
		rCopy.left.extend([
		'node [ id %d label "H." ]' % (root),
		'node [ id %d label "N" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "O" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "N." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy


		
		#For H. + O=S -> H-O-S.
		rCopy = r.clone()
		rCopy.name += "H. + O=S"
		
		rCopy.left.extend([
		'node [ id %d label "H." ]' % (root),
		'node [ id %d label "S" ]' % (offset+1),		
		'edge [ source %d target %d label "=" ]' % (offset, offset+1),
		
		])		

		rCopy.context.extend([
		'node [ id %d label "O" ]' % (offset),				
		
		])
		
		rCopy.right.extend([
		'node [ id %d label "H" ]' % (root),
		'node [ id %d label "S." ]' % (offset+1),		
		'edge [ source %d target %d label "-" ]' % (offset, offset+1),		
		'edge [ source %d target %d label "-" ]' % (root, offset),

		])
			
		yield rCopy

		
	r = RuleGen("doublebondattack")
	for r in doublemrad(r, 0, 1):
		yield r.loadRule()

doublebondattack = [a for a in doublebondattackGen()]
