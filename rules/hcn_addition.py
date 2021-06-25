# Meant to be reversible
# Modified such that it only does HCN + HCN, not RCN+HCN (less common reaction)
nitrile_addition = ruleGMLString("""rule [
    ruleID "HCN Addition to Nitriles"
    labelType "term"
    left [
        edge [ source 2 target 3 label "#" ]
        edge [ source 5 target 6 label "-" ]
    ]
    context [
        node [ id 1 label "H" ] # note that the * could even be Oxygen
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

# again, reversible
add_iminonitrile = ruleGMLString("""rule [
    ruleID "CN addition to iminonitriles"
    labelType "term"
    left [
        edge [ source 2 target 3 label "=" ]
        edge [ source 7 target 8 label "-" ]
    ]
    context [
        node [ id 1 label "*" ]
        edge [ source 1 target 2 label "-" ]
        node [ id 2 label "C" ]
        node [ id 3 label "N" ] # imine
        edge [ source 3 target 4 label "-" ]
        node [ id 4 label "H" ]
        edge [ source 2 target 5 label "-" ]
        node [ id 5 label "C" ]
        edge [ source 5 target 6 label "#" ]
        node [ id 6 label "N" ]

        node [ id 7 label "H" ]
        node [ id 8 label "C" ]
        node [ id 9 label "N" ]
        edge [ source 8 target 9 label "#" ]

    ]
    right [
        edge [ source 2 target 3 label "-" ]
        edge [ source 2 target 8 label "-" ]
        edge [ source 7 target 3 label "-" ] # convert into amine
    ]
]""")


elim_iminonitrile = add_iminonitrile.makeInverse()


add_aminodinitriles = ruleGMLString("""rule [
    ruleID "CN addition to aminodinitriles"
    left [
        edge [ source 1 target 2 label "-" ]
        edge [ source 2 target 5 label "-" ]
        edge [ source 5 target 6 label "#" ]
        edge [ source 10 target 11 label "-" ]
    ]
    context [
        node [ id 1 label "H" ]
        node [ id 2 label "C" ]
        edge [ source 2 target 3 label "-" ]
        node [ id 3 label "C" ]
        edge [ source 3 target 4 label "#" ]
        node [ id 4 label "N" ]
        
        node [ id 5 label "C" ]
        node [ id 6 label "N" ]
        edge [ source 2 target 7 label "-" ]
        node [ id 7 label "N" ]
        edge [ source 7 target 8 label "-" ]
        node [ id 8 label "*" ]
        edge [ source 7 target 9 label "-" ]
        node [ id 9 label "*" ]

        node [ id 10 label "H" ]
        node [ id 11 label "C" ]
        edge [ source 11 target 12 label "#" ]
        node [ id 12 label "N" ]
    ]
    right [
        edge [ source 2 target 5 label "=" ]
        edge [ source 5 target 6 label "-" ]
        edge [ source 1 target 6 label "-" ]
        edge [ source 10 target 6 label "-" ]
        edge [ source 5 target 11 label "-" ]
    ]
]""")

elim_aminodinitriles = add_aminodinitriles.makeInverse()

'''nitrile_addition.print()
nitrile_addition_inv.print()
add_iminonitrile.print()
elim_iminonitrile.print()
add_aminodinitriles.print()
elim_aminodinitriles.print()'''