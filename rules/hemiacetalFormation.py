# Note: this is specifically for 6 membered (including oxygen) rings; such a huge constraint
# has been put to reduce the number of products

hemiacetalFormation = [ruleGMLString("""rule [
    ruleID "Hemiacetal Formation for Rings"
    left [
        edge [ source 2 target 3 label "-" ]
        edge [ source 4 target 5 label "=" ]
    ]
    context [
        node [ id 1 label "C" ]
        node [ id 2 label "O" ]
        node [ id 3 label "H" ]
        node [ id 4 label "C" ]
        node [ id 5 label "O" ]
        edge [ source 1 target 2 label "-" ]
        # additional context to limit this to only for forming/breaking 6 membered rings
        node [ id 6 label "C" ]
        node [ id 7 label "C" ]
        node [ id 8 label "C" ]
        edge [ source 1 target 6 label "-" ]
        edge [ source 6 target 7 label "-" ]
        edge [ source 7 target 8 label "-" ]
        edge [ source 8 target 4 label "-" ]
    ]
    right [
        edge [ source 4 target 5 label "-" ]
        edge [ source 3 target 5 label "-" ]
        # close the ring by joining O and C
        edge [ source 4 target 2 label "-" ]
    ]
]
""")]