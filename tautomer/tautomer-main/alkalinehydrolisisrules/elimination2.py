betaElimination2 = [ruleGMLString("""
rule [
    ruleID "Beta Elimination2"
    left  [  
        edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 1 target 4 label "-" ]
		edge [ source 1 target 5 label "-" ]
		edge [ source 8 target 5 label "-" ]
		edge [ source 5 target 6 label "-" ]
		edge [ source 6 target 7 label "-" ]
		edge [ source 9 target 5 label "-" ]
		edge [ source 9 target 10 label "-" ]
		edge [ source 10 target 11 label "-" ]
    ]
    context [
        node [ id 1 label "C" ]
		node [ id 2 label "O" ]
		node [ id 3 label "H" ]
		node [ id 4 label "H" ]
		node [ id 5 label "C" ]
		node [ id 6 label "O" ]
		node [ id 7 label "H" ]
        node [ id 8 label "H" ]
		node [ id 9 label "C" ]
		node [ id 10 label "O" ]
        node [ id 11 label "H" ]
    ]
    right [
        edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 6 target 4 label "-" ]
		edge [ source 1 target 5 label "=" ]
		edge [ source 8 target 5 label "-" ]
		edge [ source 6 target 7 label "-" ]
		edge [ source 9 target 5 label "-" ]
		edge [ source 9 target 10 label "-" ]
		edge [ source 10 target 11 label "-" ]
	]

]""")]
print(inputRules)

postSection ("Input Rules")
for a in inputRules:
	a.print()