alpha_beta_unsaturated = ruleGMLString("""rule [
    ruleID "Alpha-Beta Unsaturated Acid Decarboxylation"
    labelType "term"
    left [
        edge [ source 4 target 6 label "-" ]
        edge [ source 6 target 8 label "-" ]
        edge [ source 8 target 9 label "-" ]
    ]
    context [
        node [ id 1 label "*" ]
        edge [ source 1 target 2 label "-" ]
        node [ id 2 label "C" ]
        edge [ source 2 target 3 label "-" ]
        node [ id 3 label "*" ]
        edge [ source 2 target 4 label "=" ]
        node [ id 4 label "C" ]
        edge [ source 4 target 5 label "-" ]
        node [ id 5 label "H" ]
        node [ id 6 label "C" ]
        edge [ source 6 target 7 label "=" ]
        node [ id 7 label "O" ]
        node [ id 8 label "O" ]
        node [ id 9 label "H" ]
    ]
    right [
        edge [ source 6 target 8 label "=" ]
        edge [ source 9 target 4 label "-" ]
    ]
]""")

beta_gamma_unsaturated = ruleGMLString("""rule [
    ruleID "Beta-gamma Unsaturated Acid Decarboxylation"
    labelType "term"
    left [
        edge [ source 6 target 9 label "-" ]
        edge [ source 9 target 11 label "-" ]
        edge [ source 11 target 12 label "-" ]
    ]
    context [
        node [ id 1 label "*" ]
        node [ id 2 label "C" ]
        node [ id 3 label "H" ]
        node [ id 4 label "C" ]
        node [ id 5 label "*" ]
        node [ id 6 label "C" ]
        node [ id 7 label "*" ]
        node [ id 8 label "*" ]
        node [ id 9 label "C" ]
        node [ id 10 label "O" ]
        node [ id 11 label "O" ]
        node [ id 12 label "H" ]
        edge [ source 1 target 2 label "-" ]
        edge [ source 2 target 3 label "-" ]
        edge [ source 2 target 4 label "=" ]
        edge [ source 4 target 5 label "-" ]
        edge [ source 4 target 6 label "-" ]
        edge [ source 6 target 7 label "-" ]
        edge [ source 6 target 8 label "-" ]
        edge [ source 9 target 10 label "=" ]
    ]
    right [
        edge [ source 9 target 11 label "=" ]
        edge [ source 12 target 6 label "-" ]
    ]
]""")