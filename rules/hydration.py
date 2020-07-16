# The purpose of this is to provide an "inverse" rule to elimination

hydration1 = [ruleGMLString("""rule [
    ruleID "Hydration"
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

# Hydration + 
hydration2 = [ruleGMLString("""rule [
    left [

    ]
    context [

    ]
    right [
        
    ]
]
""")]