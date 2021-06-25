heyns_rearrangement = [ruleGMLString("""rule [
    ruleID "Amadori/Heyns Rearrangement"
    labelType "term"
    left [
        edge [ source 1 target 4 label "-" ]
        edge [ source 5 target 7 label "=" ]
        edge [ source 6 target 8 label "-" ]
        edge [ source 6 target 10 label "-"]
        edge [ source 8 target 9 label "-" ]
    ]
    context [
        node [ id 1 label "N" ]
        node [ id 2 label "*" ] # Treat R1 as *
        node [ id 3 label "*" ] # Treat R2 as *
        node [ id 4 label "H" ]
        edge [ source 1 target 2 label "-" ]
        edge [ source 1 target 3 label "-" ]
        edge [ source 5 target 6 label "-" ]
        node [ id 5 label "C" ]
        node [ id 6 label "C" ]
        node [ id 7 label "O" ]
        node [ id 8 label "O" ]
        node [ id 9 label "H" ]
        node [ id 10 label "H"]
    ]
    right [
        edge [ source 1 target 5 label "-" ]
        edge [ source 5 target 10 label "-" ]
        edge [ source 6 target 8 label "=" ]

        # Remove H2O
        edge [ source 7 target 4 label "-" ]
        edge [ source 7 target 9 label "-" ]
    ]
    # if R1 and R2 are C then it must not have a (=O), i.e. the -NH2 shouldn't be part of an amide
    constrainAdj [ id 2 op "=" count 0 
        nodeLabels [ label "O" ]
        edgeLabels [ label "=" ]
    ]
    constrainAdj [ id 3 op "=" count 0 
        nodeLabels [ label "O" ]
        edgeLabels [ label "=" ]
    ]
    # The carbonyl should not be a part of a carboxylic acid/amide/thiocarboxylic acid
    constrainAdj [ id 5 op "=" count 0
        nodeLabels [ label "O" label "N" label "S" ]
        edgeLabels [ label "-" ]
    ]
]""")]

'''heyns_rearrangement_inverse = [ruleGMLString("""rule [
    ruleID "inverse Amadori/Heyns Rearrangement"
    left [
        edge [ source 1 target 5 label "-" ]
        edge [ source 5 target 10 label "-" ]
        edge [ source 6 target 8 label "=" ]

        # Add H2O
        edge [ source 7 target 4 label "-" ]
        edge [ source 7 target 9 label "-" ]
    ]
    context [
        node [ id 1 label "N" ]
        node [ id 2 label "*" ] # Treat R1 as *
        node [ id 3 label "*" ] # Treat R2 as *
        node [ id 4 label "H" ]
        edge [ source 1 target 2 label "-" ]
        edge [ source 1 target 3 label "-" ]
        edge [ source 5 target 6 label "-" ]
        node [ id 5 label "C" ]
        node [ id 6 label "C" ]
        node [ id 7 label "O" ]
        node [ id 8 label "O" ]
        node [ id 9 label "H" ]
        node [ id 10 label "H"]
    ]
    right [
        edge [ source 1 target 4 label "-" ]
        edge [ source 5 target 7 label "=" ]
        edge [ source 6 target 8 label "-" ]
        edge [ source 6 target 10 label "-"]
        edge [ source 8 target 9 label "-" ]
    ]
]""")]'''