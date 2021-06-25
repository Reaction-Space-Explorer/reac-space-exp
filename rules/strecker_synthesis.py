strecker_syn = ruleGMLString("""rule [
    ruleID "Strecker Synthesis"
    labelType "term"
    left [
        edge [ source 2 target 3 label "=" ]
        edge [ source 5 target 8 label "-" ]
        edge [ source 9 target 10 label "-" ]
    ]
    context [
        node [ id 1 label "*" ]
        node [ id 2 label "C" ]
        edge [ source 1 target 2 label "-" ]
        node [ id 3 label "O" ]
        edge [ source 4 target 2 label "-" ]
        node [ id 4 label "*" ]

        node [ id 5 label "N" ]
        node [ id 6 label "*" ]
        edge [ source 6 target 5 label "-" ]
        node [ id 7 label "*" ]
        edge [ source 7 target 5 label "-" ]
        node [ id 8 label "H" ]

        node [ id 9 label "H" ]
        node [ id 10 label "C" ]
        edge [ source 10 target 11 label "#" ]
        node [ id 11 label "N" ]
    ]
    right [
        edge [ source 2 target 5 label "-" ]
        edge [ source 10 target 2 label "-" ]

        #-H2O
        edge [ source 3 target 8 label "-" ]
        edge [ source 3 target 9 label "-" ]
    ]
    constrainAdj [
        id 2 op "=" count 0
        nodeLabels [ label "O" label "N" label "S" ]
        edgeLabels [ label "-" ]
    ]
    constrainAdj [
        id 5 op "=" count 0
        nodeLabels [ label "O" label "S" ]
        edgeLabels [ label "-" ]
    ]
    constrainAdj [
        id 6 op "=" count 0
        nodeLabels [ label "O" label "S" ]
        edgeLabels [ label "=" ]
    ]
    constrainAdj [
        id 7 op "=" count 0
        nodeLabels [ label "O" label "S" ]
        edgeLabels [ label "=" ]
    ]
]""")

#strecker_syn.print()