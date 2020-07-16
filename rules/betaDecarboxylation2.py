include("common.py")

def betaDecarboxylation2Gen():
	r = RuleGen("Beta-decarboxylation 2")
	r.label = "term"
	r.left.extend([
		'edge [ source 0 target 3 label "-" ]',
		'edge [ source 3 target 5 label "-" ]',
		'edge [ source 5 target 6 label "-" ]',
	])
	r.context.extend([
		'node [ id 0 label "C" ]',
		'node [ id 1 label "*" ]',
		'node [ id 2 label "*" ]',
		'edge [ source 0 target 1 label "-" ]',
		'edge [ source 0 target 2 label "-" ]',
		'node [ id 3 label "C" ]',
		'node [ id 4 label "O" ]',
		'edge [ source 3 target 4 label "=" ]',
		'node [ id 5 label "O" ]',
		'node [ id 6 label "H" ]',
		'# Variation part',
		'node [ id 10 label "C" ]',
		'edge [ source 0 target 10 label "-" ]',
	])
	r.right.extend([
		'edge [ source 3 target 5 label "=" ]',
		'edge [ source 0 target 6 label "-" ]',
	])
	# CN
	r1 = r.clone()
	r1.name += " CN"
	r1.context.extend([
		'node [ id 11 label "N" ]',
		'edge [ source 10 target 11 label "#" ]',
	])
	yield r1.loadRule()
	# CO*
	r.name += " CO"
	r.context.extend([
		'node [ id 11 label "O" ]',
		'edge [ source 10 target 11 label "=" ]',
		'node [ id 13 label "H" ]',
		'edge [ source 10 target 12 label "-" ]',
		'edge [ source 12 target 13 label "-" ]',
	])
	# COOH
	r1 = r.clone()
	r1.name += "OH"
	r1.context.extend([
		'node [ id 12 label "O" ]',
	])
	yield r1.loadRule()
	# CONH2
	r1 = r.clone()
	r1.name += "NH2"
	r1.context.extend([
		'node [ id 12 label "N" ]',
		'node [ id 14 label "H" ]',
		'edge [ source 12 target 14 label "-" ]',
	])
	yield r1.loadRule()
	# COSH
	r1 = r.clone()
	r1.name += "SH"
	r1.context.extend([
		'node [ id 12 label "S" ]',
	])
	yield r1.loadRule()
	
	

betaDecatboxylation2 = [a for a in betaDecarboxylation2Gen()]
