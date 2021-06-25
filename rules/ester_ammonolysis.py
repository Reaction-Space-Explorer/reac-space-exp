ammonolysis_ester = ruleGMLString("""rule [
    ruleID "Ammonolysis of Esters"
    labelType "term"
    left [
        edge [ source 2 target 4 label "-" ]
        edge [ source 6 target 7 label "-"]
    ]
    context [
        node [ id 1 label "*" ]
        edge [ source 1 target 2 label "-" ]
        node [ id 2 label "C" ]
        edge [ source 2 target 3 label "=" ]
        node [ id 3 label "O" ]
        node [ id 4 label "O" ]
        edge [ source 4 target 5 label "-" ]
        node [ id 5 label "C" ]

        node [ id 6 label "H" ]
        node [ id 7 label "N" ]
        edge [ source 7 target 8 label "-" ]
        node [ id 8 label "*" ]
        edge [ source 7 target 9 label "-" ]
        node [ id 9 label "*" ]
    ]
    right [
        edge [ source 4 target 6 label "-" ]
        edge [ source 7 target 2 label "-" ]
    ]
]""")

#ammonolysis_ester.print()