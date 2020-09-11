# This is meant to be irreversible
nitrile_hydrol = ruleGMLString("""rule [
    ruleID "Nitrile Hydrolysis"
    labelType "term"
    left [
        edge [ source 2 target 3 label "#" ]
        edge [ source 4 target 5 label "-" ]
        edge [ source 4 target 6 label "-" ]
    ]
    context [
        node [ id 1 label "*" ]
        node [ id 2 label "C" ]
        edge [ source 1 target 2 label "-" ]
        node [ id 3 label "N" ]

        node [ id 4 label "O" ]
        node [ id 5 label "H" ]
        node [ id 6 label "H" ]
    ]
    right [
        edge [ source 2 target 4 label "=" ]
        edge [ source 2 target 3 label "-" ]
        edge [ source 3 target 5 label "-" ]
        edge [ source 3 target 6 label "-" ]
    ]
]
""")