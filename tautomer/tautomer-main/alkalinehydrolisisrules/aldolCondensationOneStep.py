AldolCondensationOneStep = [ruleGMLString("""
rule [
    ruleID "Aldol Condensation One Step"
    left  [  
        edge [ source 2 target 3 label "-" ]
		edge [ source 2 target 5 label "-" ]
		edge [ source 1 target 3 label "-" ]
		edge [ source 3 target 4 label "-" ]
		edge [ source 4 target 6 label "-" ]
		edge [ source 7 target 9 label "-" ]
		edge [ source 8 target 9 label "=" ]
		edge [ source 10 target 9 label "-" ]
    ]
    context [
        node [ id 1 label "H" ]
		node [ id 2 label "O" ]
		node [ id 3 label "C" ]
		node [ id 4 label "C" ]
		node [ id 5 label "H" ]
		node [ id 6 label "H" ]
		node [ id 7 label "H" ]
		node [ id 8 label "O" ]
		node [ id 9 label "C" ]
		node [ id 10 label "C" ]
    ]
    right [
        edge [ source 2 target 3 label "=" ]
		edge [ source 1 target 3 label "-" ]
		edge [ source 9 target 7 label "-" ]
		edge [ source 3 target 4 label "-" ]
		edge [ source 8 target 5 label "-" ]
		edge [ source 9 target 4 label "-" ]
		edge [ source 6 target 8 label "-" ]
		edge [ source 10 target 9 label "-" ]
	]

]""")]
print(inputRules)

postSection ("Input Rules")
for a in inputRules:
	a.print()