include("common.py")
#A realizar una mejor busqueda de terminaciones, ver tambien la auto reaccion entre radicales 
def hydroxylalcaneGen():
	def ohalcane(r, root, offset):
		
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

				
	r = RuleGen("OHRadAlcanes")
	for r in ohalcane(r, 0, 1):
		yield r.loadRule()

ohradalcanes = [a for a in hydroxylalcaneGen()]
