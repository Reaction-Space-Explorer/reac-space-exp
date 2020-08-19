include("common.py")

# Some rules have a missing adjConstraint due to which they can't be inverted
#config.rule.ignoreConstraintsDuringInversion = True

if True:
	# The ones that were created by Rana and Aayush
	include("aldolCondensation.py")
	include("benzilicAcidRearrangement.py")
	include("elimination1.py")
	include("elimination2.py")
	include("retroAldol.py")
	include("ketoEnolisation.py")
	include("hemiacetalFormation.py")
	include("hydration.py")
	include("heyns_rearrangement.py")
	# The following were taken from the HCN folder
	include("esterFormationHydrolysisExchange.py")
	include("cannizarro1.py")
	include("cannizarro2.py")
	include("alkeneAdditionElimination.py")
	include("alkyneAddition.py")
	include("knoevenagel.py")
	include("damnScission.py")
	include("cyanogenAmmonolysis.py")
	include("cyanamideHydration.py")
	include("carbamylation.py")
	include("cyanamidation.py")
	include("mannich.py")
	include("amideFormationHydrolysis.py")
	include("cyanateFromCyanogen.py")
	include("betaDecarboxylation2.py")
	include("transamination.py")
	include("alphaKetoAcidDecarboxylation.py")
	include("ketoneCatalyzedDecarboxylation.py")
	include("amadoriRearrangement.py")
	include("streckerDegradation.py")
	include("amideExchange.py")
	include("thialThiene.py")
	include("sulfideToNitrileAddition.py")
	include("amidineToAmideHydrolysis.py")
	include("thioamidineToAmideHydrolysis.py")
	include("ketoneAldehydeToThioketoneThial.py")
	include("ringClosure.py")
	include("hcnToKetone.py")
	include("imineToCarbonyl.py")
	include("imineToEnamine.py")
	include("streckerDegradationDicarbonyl.py")
	include("schiffTautomerization.py")
	include("michaelAddition.py")

canNotBeInvertedYet = []
loaded = list(inputRules)
for a in loaded:
	lt = a.labelType
	if lt is None:
		def checkLabel(l):
			if l == "*":
				print("Missing labelType in rule:", a.name)
				a.print()
				assert False
		for v in a.vertices:
			vl = v.left
			vr = v.right
			if not vl.isNull():
				checkLabel(vl.stringLabel)
			if not vr.isNull():
				checkLabel(vr.stringLabel)
		for e in a.edges:
			el = e.left
			er = e.right
			if not el.isNull():
				checkLabel(el.stringLabel)
			if not er.isNull():
				checkLabel(er.stringLabel)
	try:
		inv = a.makeInverse()
		if inv.isomorphism(a) == 0:
			inputRules.append(inv)
	except LogicError as e:
		print(a.name, "can not be inverted yet.")
		canNotBeInvertedYet.append(a)
