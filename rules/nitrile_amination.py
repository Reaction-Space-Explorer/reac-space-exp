# amination of RCN

nitr_amination = ruleGMLString("""rule [
    ruleID "Nitrile Amination"
    labelType "term"
    left [
        edge [ source 2 target 3 label "#" ]
        edge [ source 5 target 6 label "-" ]
    ]
    context [
        node [ id 1 label "*" ]
        edge [ source 1 target 2 label "-" ]
        node [ id 2 label "C" ]
        node [ id 3 label "N" ]

        node [ id 4 label "*" ]
        edge [ source 4 target 5 label "-" ]
        node [ id 5 label "N" ]
        node [ id 6 label "H" ]
        edge [ source 5 target 7 label "-" ]
        node [ id 7 label "*" ]
    ]
    right [
        edge [ source 2 target 3 label "=" ] # convert into imine
        edge [ source 5 target 2 label "-" ] # C-N
        edge [ source 6 target 3 label "-" ] # N-H
    ]
]""")