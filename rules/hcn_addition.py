# Meant to be reversible
nitrile_addition = ruleGMLString("""rule [
    ruleID "HCN Addition to Nitriles"
    labelType "term"
    left [
        edge [ source 2 target 3 label "#" ]
        edge [ source 5 target 6 label "-" ]
    ]
    context [
        node [ id 1 label "*" ] # note that the * could even be Oxygen
        edge [ source 1 target 2 label "-" ]
        node [ id 2 label "C" ]
        node [ id 3 label "N" ]

        node [ id 4 label "N" ]
        edge [ source 4 target 5 label "#" ]
        node [ id 5 label "C" ]
        node [ id 6 label "H"]
    ]
    right [
        edge [ source 2 target 5 label "-" ]
        edge [ source 3 target 6 label "-" ]
        edge [ source 2 target 3 label "=" ]
    ]
]""")

nitrile_addition_inv = nitrile_addition.makeInverse()