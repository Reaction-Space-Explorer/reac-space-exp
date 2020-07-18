# The purpose of this is to provide an "inverse" rule to elimination

# Hydration + 
hydration1 = [ruleGMLString("""rule [
    ruleID "Hydration of C(=O)C"
    left [
        edge [ source 1 target 2 label "=" ]
        edge [ source 3 target 4 label "-" ]
        edge [ source 5 target 6 label "-" ]
    ]
    context [
        node [ id 1 label "C" ]
        node [ id 2 label "O" ]
        node [ id 3 label "C" ]
        node [ id 4 label "H" ]
        #H2O
        node [ id 5 label "O" ]
        node [ id 6 label "H" ]
        node [ id 7 label "H" ]

        edge [ source 1 target 3 label "-" ]
        edge [ source 5 target 7 label "-" ]
    ]
    right [
        edge [ source 1 target 2 label "-" ]
        edge [ source 1 target 4 label "-" ]
        edge [ source 2 target 6 label "-" ]
        edge [ source 5 target 3 label "-" ]
    ]
    constrainAdj [ id 1 op "=" count 0 
        edgeLabels [ label "-" ]
        nodeLabels [ label "N" label "S" label "O" ]
    ]
    constrainAdj [ id 3 op "=" count 0
        edgeLabels [ label "=" label "-" ] # avoid gem diols as well
        nodeLabels [ label "O" label "N" label "S" ]
    ]
]
""")]

hydration2 = [ruleGMLString("""rule [
    ruleID "Hydration of C=C(O)"
    left [
        edge [ source 1 target 2 label "=" ]
        edge [ source 3 target 4 label "-" ]
    ]
    context [
        node [ id 1 label "C" ]
        node [ id 2 label "C" ]

        node [ id 3 label "H" ]
        node [ id 4 label "O" ]
        node [ id 5 label "H" ]
        edge [ source 4 target 5 label "-" ]
    ]
    right [
        edge [ source 1 target 2 label "-" ]
        edge [ source 4 target 2 label "-" ] # join -OH
        edge [ source 3 target 1 label "-" ] # join -H
    ]
    constrainAdj [ id 2 op "=" count 0
        nodeLabels [ label "N" label "S" ]
        edgeLabels [ label "-" label "=" ]
    ]
    # It should not form carboxylic acids (avoid adding -OH to carbons with =O )
    constrainAdj [ id 2 op "=" count 0
        nodeLabels [ label "O" ]
        edgeLabels [ label "=" ]
    ]
]
""")]