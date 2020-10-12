cyanohydrin_reac = ruleGMLString("""rule [
    ruleID "Cyanohydrin Formation"
    labelType "term"
    left [
        edge [ source 5 target 6 label "-" ]
        edge [ source 2 target 3 label "=" ]
    ]
    context [
        node [ id 1 label "*" ]
        edge [ source 1 target 2 label "-" ]
        node [ id 2 label "C" ]
        node [ id 3 label "O" ]
        edge [ source 2 target 4 label "-" ]
        node [ id 4 label "*" ]

        node [ id 5 label "H" ]
        node [ id 6 label "C" ]
        edge [ source 6 target 7 label "#" ]
        node [ id 7 label "N" ]
    ]
    right [
        edge [ source 2 target 3 label "-" ]
        edge [ source 3 target 5 label "-" ]
        edge [ source 2 target 6 label "-" ]
    ]
    constrainAdj [
        id 2 op "=" count 0
        nodeLabels [ label "O" label "N" label "S" ]
        edgeLabels [ label "-" ]
    ]
    
]""")
'''constrainAdj [ 
        id 1 op "=" count 0
        nodeLabels [ label "O" label "N" label "S" ]
        edgeLabels [ label "=" ]
    ]'''
