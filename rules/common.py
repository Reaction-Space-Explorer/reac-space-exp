class RuleGen(object):
	def __init__(self, name):
		self.name = name
		self.label = None
		self.left = []
		self.context = []
		self.right = []
		self.constraints = []

	def clone(self):
		r = RuleGen(self.name)
		r.label = self.label
		r.left = list(self.left)
		r.context = list(self.context)
		r.right = list(self.right)
		r.constraints = list(self.constraints)
		return r

	def loadRule(self):
		def addSection(g, name, lines):
			g += "\t%s [\n" % name
			for l in lines:
				g += "\t\t%s\n" % l
			return g + "\t]\n"
		gml = 'rule [\n\truleID "%s"\n' % self.name
		if self.label:
			gml += "labelType \"%s\"\n" % self.label
		gml = addSection(gml, "left", self.left)
		gml = addSection(gml, "context", self.context)
		gml = addSection(gml, "right", self.right)
		for l in self.constraints:
			gml += "\t%s\n" % l
		gml += "]\n"
		return ruleGMLString(gml)


def checkRuleRelation(mono=True, doPrint=False):
	for i in range(len(inputRules)):
		for j in range(i + 1, len(inputRules)):
			a = inputRules[i]
			b = inputRules[j]
			if a.isomorphism(b):
				print("Isomorphism:\n\t%s\n\t%s\n" % (a.name, b.name))
				if doPrint:
					postSection("Isomorphism")
					a.print()
					b.print()
				continue
			if mono:
				if a.monomorphism(b):
					print("Monomorphism:\n\t%s\n\t%s\n" % (a.name, b.name))
					if doPrint:
						postSection("Monomorphism")
						a.print()
						b.print()
					continue
				if b.monomorphism(a):
					print("Monomorphism:\n\t%s\n\t%s\n" % (b.name, a.name))
					if doPrint:
						postSection("Monomorphism")
						b.print()
						a.print()
					continue


def attach_H_C(r, root, offset, withCarbonylConstraint=False):
	r.context.extend([
		"# H or C on %d" % root,
		'edge [ source %d target %d label "-" ]' % (root, offset),
	])
	for s in ["H", "C"]:
		rCopy = r.clone()
		rCopy.name += ", " + s
		rCopy.context.append('node [ id %d label "%s" ]' % (offset, s))
		if s == 'C':
			rCopy.constraints.extend([
				'constrainAdj [',
				'	id %d op "=" count 0' % offset,
				'	nodeLabels [ label "O" ]',
				'	edgeLabels [ label "=" ]',
				']',
			])
		yield rCopy


def attach_2_H_C(r, root, offset, withCarbonylConstraint=False):
	r.context.extend([
		"# 2, H or C on %d" % root,
	])
	for s1, s2 in [("H", "H"), ("H", "C"), ("C", "C")]:
		rCopy = r.clone()
		rCopy.name += ", " + s1 + s2
		rCopy.context.extend([
			'edge [ source %d target %d label "-" ]' % (root, offset),
			'node [ id %d label "%s" ]' % (offset, s1),
			'edge [ source %d target %d label "-" ]' % (root, offset + 1),
			'node [ id %d label "%s" ]' % (offset + 1, s2),
		])
		if s1 == 'C':
			rCopy.constraints.extend([
				'constrainAdj [',
				'	id %d op "=" count 0' % offset,
				'	nodeLabels [ label "O" ]',
				'	edgeLabels [ label "=" ]',
				']',
			])
		if s2 == 'C':
			rCopy.constraints.extend([
				'constrainAdj [',
				'	id %d op "=" count 0' % (offset + 1),
				'	nodeLabels [ label "O" ]',
				'	edgeLabels [ label "=" ]',
				']',
			])
		yield rCopy


def attach_EWG(r, root, offset):
	r.context.extend([
		"# EWG on %d" % root,
	])
	
	rCopy = r.clone()
	rCopy.name += ", CN"
	rCopy.context.extend([
		'edge [ source %d target %d label "-" ]' % (root, offset),
		'node [ id %d label "C" ]' % offset,
		'node [ id %d label "N" ]' % (offset + 1),
		'edge [ source %d target %d label "#" ]' % (offset, offset + 1),
	])
	yield rCopy

	rCopy = r.clone()
	rCopy.name += ", COOR"
	rCopy.context.extend([
		'edge [ source %d target %d label "-" ]' % (root, offset),
		'node [ id %d label "C" ]' % offset,
		'node [ id %d label "O" ]' % (offset + 1),
		'edge [ source %d target %d label "=" ]' % (offset, offset + 1),
		'node [ id %d label "O" ]' % (offset + 2),
		'edge [ source %d target %d label "-" ]' % (offset, offset + 2),
	])
	for rCopy in attach_H_C(rCopy, offset + 2, offset + 3):
		yield rCopy

	rCopy = r.clone()
	rCopy.name += ", CONRR"
	rCopy.context.extend([
		'edge [ source %d target %d label "-" ]' % (root, offset),
		'node [ id %d label "C" ]' % offset,
		'node [ id %d label "O" ]' % (offset + 1),
		'edge [ source %d target %d label "=" ]' % (offset, offset + 1),
		'node [ id %d label "N" ]' % (offset + 2),
		'edge [ source %d target %d label "-" ]' % (offset, offset + 2),
	])
	for rCopy in attach_2_H_C(rCopy, offset + 2, offset + 3):
		yield rCopy

	rCopy = r.clone()
	rCopy.name += ", NRR"
	rCopy.context.extend([
		'edge [ source %d target %d label "-" ]' % (root, offset),
		'node [ id %d label "N" ]' % offset,
	])
	for rCopy in attach_2_H_C(rCopy, offset, offset + 1):
		yield rCopy

	rCopy = r.clone()
	rCopy.name += ", OR"
	rCopy.context.extend([
		'edge [ source %d target %d label "-" ]' % (root, offset),
		'node [ id %d label "O" ]' % offset,
	])
	for rCopy in attach_H_C(rCopy, offset, offset + 1):
		yield rCopy

	rCopy = r.clone()
	rCopy.name += ", C(=O)R"
	rCopy.context.extend([
		'edge [ source %d target %d label "-" ]' % (root, offset),
		'node [ id %d label "C" ]' % offset,
		'node [ id %d label "O" ]' % (offset + 1),
		'edge [ source %d target %d label "=" ]' % (offset, offset + 1),
	])
	for rCopy in attach_H_C(rCopy, offset, offset + 2):
		yield rCopy

def attach_Nu(r, offset, allow_o=True):
	'''
	allow_o -- whether or not to consider O for attachment or not. I added this as a
	quick fix to avoid attaching O in the alkene-addition elimination rule (which was
	producing enols). It is passed as False when that rule is generated.
	'''
	rOrig = r
	for k in "ONSC":
		r = rOrig.clone()
		r.name += ", %s" % k
		r.context.extend([
			'node [ id %d label "%s" ]' % (offset, k),
		])
		if k == 'O':
			if allow_o == True:
				for r in attach_H_C(r, offset, offset + 1):
					yield r
		elif k == 'S':
			for r in attach_H_C(r, offset, offset + 1):
				yield r
		elif k == 'N':
			r.context.extend([
				'# TODO: constrain: N not part of an amide'
			])
			for r in attach_2_H_C(r, offset, offset + 1):
				yield r
		elif k == 'C':
			r.context.extend([
				'node [ id %d label "N" ]' % (offset + 1),
				'edge [ source %d target %d label "#" ]' % (offset, offset + 1),
			])
			yield r
		else:
			assert False









