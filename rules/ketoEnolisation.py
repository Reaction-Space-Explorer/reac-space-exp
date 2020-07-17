ketoEnolisation = [ruleGMLString("""rule [
    ruleID "Keto-enol migration twice" # to avoid enols (and tautomers in general)
    left[
        edge [ source 1 target 2 label "=" ]
        edge [ source 3 target 4 label "-" ]
    ]
    context [
        node [ id 1 label "C" ]
        node [ id 2 label "O" ]
        node [ id 3 label "C" ]
        node [ id 4 label "O" ]
        node [ id 5 label "H" ]
        edge [ source 4 target 5 label "-" ]
        edge [ source 1 target 3 label "-" ]
    ]
    right [
        edge [ source 1 target 4 label "-" ]
        edge [ source 3 target 2 label "=" ]
    ]
]
""")]