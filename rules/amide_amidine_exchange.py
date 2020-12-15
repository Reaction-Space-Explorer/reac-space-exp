"""
So far, these are only intramolecular amide/amidine exchange.
"""

amide_amidine_ex_5 = ruleGMLString("""rule [
    ruleID "Intramolecular Amide-Amidine Exchange, 5"
    labelType "term"
    left [
        edge [ source 2 target 3 label "-" ]
        edge [ source 6 target 8 label "-" ]
    ]
    context [
        node [ id 1 label "*" ] # supposed to be C/H
        node [ id 2 label "N" ]
        edge [ source 1 target 2 label "-" ]
        node [ id 3 label "H" ]
        node [ id 4 label "C" ]
        edge [ source 2 target 4 label "-" ]
        node [ id 5 label "C" ]
        edge [ source 4 target 5 label "*" ]
        node [ id 6 label "N" ]
        edge [ source 5 target 6 label "-" ]
        node [ id 7 label "*" ] # C/H
        edge [ source 6 target 7 label "-" ]
        node [ id 8 label "C" ]
        node [ id 9 label "*" ] # =O or =N
        edge [ source 8 target 9 label "=" ]
        node [ id 10 label "*" ] # C/N/H
        edge [ source 8 target 10 label "-" ]
    ]
    right [
        edge [ source 3 target 6 label "-" ]
        edge [ source 2 target 8 label "-" ]
    ]
    # make sure 1 isn't bonded to =N or =O (shouldn't be part of an amide/amidine)
    constrainAdj [ id 1 op "=" count 0 
        nodeLabels [ label "O" label "N" ]
        edgeLabels [ label "-" ]
    ]
    # make sure 6 isn't bonded to -O or -N, to make sure 7 isnt anything other than C/H
    constrainAdj [ id 6 op "=" count 0 
        edgeLabels [ label "-" ]
        nodeLabels [ label "O" label "N" ]
    ]
    # make sure 8 is bonded to one =O or =N, to ensure 9 is =O or =N
    constrainAdj [ id 8 op "=" count 1 
        edgeLabels [ label "=" ]
        nodeLabels [ label "O" label "N" ]
    ]
    # make sure 8 isn't bonded to -O, as 10 shouldn't be anything other than C/N/H
    constrainAdj [ id 8 op "=" count 0
        edgeLabels [ label "-" ]
        nodeLabels [ label "O" ]
    ]
]""")

p = GraphPrinter()
p.withColour = True
p.withIndex = True

amide_amidine_ex_5.print(p)

amide_amidine_ex_6 = ruleGMLString("""rule [
    ruleID "Intramolecular Amide-Amidine Exchange, 6"
    labelType "term"
    left [
        edge [ source 2 target 3 label "-" ]
        edge [ source 7 target 9 label "-" ]
    ]
    context [
        node [ id 1 label "*" ] # supposed to be C/H
        node [ id 2 label "N" ]
        edge [ source 1 target 2 label "-" ]
        node [ id 3 label "H" ]
        node [ id 4 label "C" ]
        edge [ source 2 target 4 label "-" ]
        node [ id 5 label "*" ]
        edge [ source 4 target 5 label "*" ]
        node [ id 6 label "*" ]
        edge [ source 5 target 6 label "*" ]
        node [ id 7 label "N" ]
        edge [ source 6 target 7 label "-" ]
        node [ id 8 label "*" ] # C/H
        edge [ source 7 target 8 label "-" ]
        node [ id 9 label "C" ]
        node [ id 10 label "*" ] # =O or =N
        edge [ source 9 target 10 label "=" ]
        node [ id 11 label "*" ] # C/N/H
        edge [ source 9 target 11 label "-" ]
    ]
    right [
        edge [ source 3 target 7 label "-" ]
        edge [ source 2 target 9 label "-" ]
    ]
    # make sure 1 isn't bonded to =N or =O (shouldn't be part of an amide/amidine)
    constrainAdj [ id 1 op "=" count 0 
        nodeLabels [ label "O" label "N" ]
        edgeLabels [ label "-" ]
    ]
    # make sure 7 isn't bonded to -O or -N, to make sure 8 isnt anything other than C/H
    constrainAdj [ id 7 op "=" count 0 
        edgeLabels [ label "-" ]
        nodeLabels [ label "O" label "N" ]
    ]
    # make sure 9 is bonded to one =O or =N, to ensure 10 is =O or =N
    constrainAdj [ id 9 op "=" count 1 
        edgeLabels [ label "=" ]
        nodeLabels [ label "O" label "N" ]
    ]
    # make sure 9 isn't bonded to -O, as 11 shouldn't be anything other than C/N/H
    constrainAdj [ id 9 op "=" count 0
        edgeLabels [ label "-" ]
        nodeLabels [ label "O" ]
    ]
]""")

amide_amidine_ex_6.print(p)