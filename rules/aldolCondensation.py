aldolCondensation = [ruleGMLString("""rule [
    ruleID "Aldol Condensation"
    left [
        edge [ source 3 target 4 label "-" ]
        edge [ source 6 target 7 label "=" ]
    ]
    context [
        node [ id 1 label "C" ]
        node [ id 2 label "O" ]
        node [ id 3 label "C" ]
        node [ id 4 label "H" ]
        node [ id 5 label "H" ] # there need to be at least 2 alpha-H on one of the carbonyls
        edge [ source 1 target 2 label "=" ]
        edge [ source 1 target 3 label "-" ]
        edge [ source 5 target 3 label "-" ]
        node [ id 6 label "C" ]
        node [ id 7 label "O" ]
    ]
    right [
        edge [ source 4 target 7 label "-" ] # convert O into OH
        edge [ source 3 target 6 label "-" ] # join the carbons
        edge [ source 7 target 6 label "-" ] # join the OH to one of the carbons
    ]
    # The carbonyl shouldn't be a part of an amide/thioamide/carboxylic acid
    constrainAdj [ id 1 op "=" count 0
        nodeLabels [ label "N" label "S" label "O" ]
        edgeLabels [ label "-" ]
    ]
    constrainAdj [ id 6 op "=" count 0
        nodeLabels [ label "N" label "S" label "O" ]
        edgeLabels [ label "-" ]
    ]
]
""")]
# The rule below does beta elimination in the same step, try to not do that
'''aldolCondensation = [ruleGMLString("""rule [
    ruleID "Aldol Condensation"
    left [
        edge [ source 1 target 2 label "=" ]
        edge [ source 5 target 6 label "-" ]
        edge [ source 5 target 7 label "-" ]
    ]
    context [
        node [ id 1 label "C" ]
        node [ id 2 label "O" ]

        # a different carbonyl group
        node [ id 3 label "C" ]
        node [ id 4 label "O" ]
        edge [ source 3 target 4 label "=" ]
        node [ id 5 label "C" ] # at least one of the carbonyls must have 2 alpha-hydrogens
        node [ id 6 label "H" ]
        node [ id 7 label "H" ]
        edge [ source 3 target 5 label "-" ]
    ]
    right [
        # H2O
        edge [ source 2 target 6 label "-" ]
        edge [ source 2 target 7 label "-" ]

        edge [ source 1 target 5 label "=" ]
    ]
]
""")]'''