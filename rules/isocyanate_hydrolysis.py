isocyanate_hydro = ruleGMLString("""rule [
    ruleID "Isocyanate Hydrolysis"
    labelType "term"
    left [
        edge [ source 2 target 3 label "=" ]
        edge [ source 6 target 5 label "-" ]
        edge [ source 6 target 7 label "-" ]
    ]
    context [
        node [ id 1 label "H" ]
        edge [ source 1 target 2 label "-" ]
        node [ id 2 label "N" ]
        node [ id 3 label "C" ]
        edge [ source 3 target 4 label "=" ]
        node [ id 4 label "O" ]

        node [ id 5 label "H" ]
        node [ id 6 label "O" ]
        node [ id 7 label "H" ]
    ]
    right [
        edge [ source 2 target 5 label "-" ]
        edge [ source 2 target 7 label "-" ]
        edge [ source 3 target 6 label "=" ]
    ]
]""")

#isocyanate_hydro.print()