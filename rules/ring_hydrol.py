include("common.py")

def create_ring_hydro(ring_size):
    r = RuleGen(f"{ring_size} Membered Ring Exoamine Hydrolysis")
    r.label = "term"
    r.left.extend([
        'edge [ source 1 target 2 label "-" ]',
        'edge [ source 100 target 101 label "-" ]',
    ])
    r.context.extend([
        'node [ id 1 label "C" ]',
        'node [ id 2 label "N" ]',
        'node [ id 3 label "*" ]',
        'edge [ source 2 target 3 label "-" ]',
        'node [ id 4 label "*" ]',
        'edge [ source 2 target 4 label "-" ]',
        'node [ id 5 label "N" ]',
        'edge [ source 1 target 5 label "=" ]',
        'node [ id 6 label "*" ]',
        'edge [ source 5 target 6 label "*" ]',
        'node [ id 7 label "*" ]',
        'edge [ source 6 target 7 label "*" ]',
        'node [ id 8 label "*" ]',
        'edge [ source 7 target 8 label "*" ]',
        'node [ id 100 label "H" ]',
        'node [ id 101 label "O" ]',
        'node [ id 102 label "H" ]',
        'edge [ source 101 target 102 label "-" ]'
    ])
    r.right.extend([
        'edge [ source 1 target 101 label "-" ]',
        'edge [ source 100 target 2 label "-" ]'
    ])
    if ring_size == 5:
        r.context.extend(['edge [ source 8 target 1 label "*" ]'])
    elif ring_size == 6:
        r.context.extend([
            'node [ id 9 label "*" ]',
            'edge [ source 8 target 9 label "*" ]',
            'edge [ source 9 target 1 label "*" ]'
        ])
    return r.loadRule()

ring_hydro_rules = [create_ring_hydro(size) for size in (5,6)]

for r in ring_hydro_rules:
    r.print()