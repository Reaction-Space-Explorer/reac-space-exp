include("common.py")

def ringClosureGen():
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
ringClosure = [a for a in ringClosureGen()]
