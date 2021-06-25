forbidden = [graphGMLString("""graph [
   node [ id 0 label "O" ]
   node [ id 1 label "O" ]
   node [ id 2 label "O" ]
   edge [ source 0 target 1 label "*" ]
   edge [ source 1 target 2 label "*" ]
   ]
   """), graphGMLString("""graph [
   node [ id 0 label "O" ]
   node [ id 1 label "O" ]
   node [ id 2 label "O." ]
   edge [ source 0 target 1 label "*" ]
   edge [ source 1 target 2 label "*" ]
   ]
   """),graphGMLString("""graph [
   node [ id 0 label "O" ]
   node [ id 1 label "O." ]
   node [ id 2 label "O" ]
   edge [ source 0 target 1 label "*" ]
   edge [ source 1 target 2 label "*" ]
   ]
   """)]

