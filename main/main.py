import os
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from molvs import tautomer

include(os.path.abspath(os.path.join('..', 'rules/all.py')))
#postChapter('Alkaline Glucose Degradation')

# starting molecule
#some_molecule = smiles('NCC=CO', name = 'starting molecule')
#some_molecule.print()
glucose = smiles('OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O', name='Glucose')

# Forbidden substructures
# Three and four membeered rings are unstable
oxirane = graphGMLString("""graph [
   node [ id 0 label "C" ]
   node [ id 1 label "C" ]
   node [ id 2 label "O" ]
   edge [ source 0 target 1 label "-"]
   edge [ source 0 target 2 label "-"]
   edge [ source 1 target 2 label "-"]
]""", name='Oxirane')

aziridine = graphGMLString("""graph [
   node [ id 0 label "C" ]
   node [ id 1 label "C" ]
   node [ id 2 label "N" ]
   edge [ source 0 target 1 label "-"]
   edge [ source 0 target 2 label "-"]
   edge [ source 1 target 2 label "-"]
]""", name='aziridine')
forbidden = [smiles('C1CC1', name='cyclopropane'), smiles('C1CCC1', name = 'cyclobutane'),
            oxirane, aziridine]

# make sure these don't get passed as an input
for fb in forbidden:
    inputGraphs.remove(fb)
strat = (
    addSubset(inputGraphs)
    >> rightPredicate[
        lambda derivation: all (g.exactMass <= 500
         and (g.monomorphism(fb) == 0 for fb in forbidden) for g in derivation.right)
    ] (inputRules)#repeat[1]()
)

# molVs smiles might not be identical to rdkit so let's use RDKit instead of molVs
enum = rdMolStandardize.TautomerEnumerator()
#enum = tautomer.TautomerEnumerator()

def clean_taut(dg, dg_execute):
    subset = dg_execute.subset
    universe = dg_execute.universe
    # The following line may seem redundant but it really isn't
    graphs = [g for g in subset]
    mod_subset_smiles = [g.smiles for g in subset]
    mod_dgverts = {v.graph: v for v in dg.vertices}
    rdkit_subset_smiles = []
    #vertex_dg_vertex = {v.graph: v for v in dg.vertices if v.graph in subset}

    taut_tautset_dict = {}
    to_remove = set()
    # Populate all possible tautomers and add to the above empty dict
    for mod_smiles in mod_subset_smiles:
        mol = Chem.MolFromSmiles(mod_smiles)
        # convert into rdkit canonical smiles
        rdkit_smiles = Chem.MolToSmiles(mol)
        rdkit_subset_smiles.append(rdkit_smiles)
        # generate Mol objects of all possible tautomers
        all_tauts = enum.Enumerate(mol)
        # convert into smiles
        all_taut_smiles = tuple(Chem.MolToSmiles(taut) for taut in all_tauts)
        #print(all_taut_smiles)
        taut_tautset_dict[mod_smiles] = all_taut_smiles

    list_tautset = tuple(taut_tautset_dict.values())
    class_ids = {}
    # create initial class ids
    for i in range(len(list_tautset)):
        class_ids.update({list_tautset[i]: i})

    # Make common taut classes
    # This is worse than O(n^2) but let's see how it goes
    for i in range(len(list_tautset)):
        for j in range(i+1, len(list_tautset)):
            for item in list_tautset[i]:
                if item in list_tautset[j]:
                    #print(f'Classes {i} and {j} have common elements')
                    #class_ids[list_tautset[i]] = class_ids[list_tautset[i]]
                    class_ids[list_tautset[j]] = class_ids[list_tautset[i]]
                    #print('Updated class ids to ', class_ids[list_tautset[i]])
                    break
    # Now, since we're done finding identical class ids
    dict_observed_tauts = {}
    for taut_class, class_id in class_ids.items():
        for i in range(len(rdkit_subset_smiles)):
            if rdkit_subset_smiles[i] in taut_class:
                # the smiles of the mod graph
                if class_id in dict_observed_tauts.keys():
                    dict_observed_tauts[class_id].append(mod_subset_smiles[i])
                else :
                    dict_observed_tauts[class_id] = [mod_subset_smiles[i]]
    print('{0} unique tautomer classes in {1} molecules'.format(len(dict_observed_tauts.keys())
                                , len(subset)))
    for taut_class in dict_observed_tauts.values():
        if len(taut_class) > 1:
            for i in range(1, len(taut_class)):
                to_remove.add(taut_class[i])
    for item in to_remove:
        subset.remove(graphs[mod_subset_smiles.index(item)])
        universe.remove(graphs[mod_subset_smiles.index(item)])          
    #for i in range(len(rdkit_subset_smiles)):

    # Merge tautomer classes
    # Now score each tautomer in each class

    return subset, universe
# What if none of the existing tautomers is the most stable one?
# We keep the most stable out of the ones in the network since we can't modify graphs in MOD
# enumerate all possible tautomers for each molecule
#   find if tautomeric pairs exist for each of them
#       if true: sort the tautomers and remove the ones with the lowest score
# items are to be removed from both subset and universe
# return the cleaned tautomers once done.

# The code became really ugly while optimizing for performance using dictionaries
'''def clean_taut(dg, dg_execute):
    subset = dg_execute.subset
    universe = dg_execute.universe

    dg_vertices_dict = {v.graph: v for v in dg.vertices if v.graph in subset}
    # A dictionary of Graph objects and the corresponding SMILES string produced by MOD
    subset_mod_smiles = {g: g.smiles for g in subset}
    inv_smiles_graph = dict(zip(subset_mod_smiles.values(), subset_mod_smiles.keys()))
    # dictionary value will later be updated to its canonical smiles
    mod_rdkit_smiles = {g.smiles: None for g in subset}
    # This dictionary will contain {canonical rdkit smiles: [list of taut smiles]}
    smiles_tauts_dict = {}
    # dictionary of canonical smiles of the isomers that need to be removed mapped to
    # the DGVertex of the one that needs to be kept
    to_remove = {}
    # for each molecule in the subset
    for mod_smiles in subset_mod_smiles.values():
        original_mol = Chem.MolFromSmiles(mod_smiles)
        # Update the dictionary value to an RDKit canonical smiles for later comparison
        mod_rdkit_smiles[mod_smiles] = Chem.MolToSmiles(original_mol)
        # enumerate all tautomers for this molecule
        all_tauts = enum.Enumerate(original_mol)
        # convert into smiles
        all_smiles = tuple(Chem.MolToSmiles(taut) for taut in all_tauts)
        smiles_tauts_dict.update({mod_rdkit_smiles[mod_smiles]: all_smiles})
    rdkit_mod_smiles = dict(zip(mod_rdkit_smiles.values(), mod_rdkit_smiles.keys()))

    # A dictionary containing {the tuple containing all possible tautomers:
    #  tautomers that were seen in the network} for each possible tautomer class
    taut_pair_dict = {}
    # Fill this dictionary
    for smiles, taut_list in smiles_tauts_dict.items():
        if taut_list in list(taut_pair_dict.keys()):
            taut_pair_dict[taut_list].append(smiles)
        else:
            taut_pair_dict.update({taut_list: [smiles]})
    # Debugging this region
    f = open('debug.txt', 'w')
    for i in range(len(taut_pair_dict.keys())-1):
        for j in range(1,len(taut_pair_dict.keys())):
            for item in list(taut_pair_dict.keys())[i]:
                if item in list(taut_pair_dict.keys())[j] and list(taut_pair_dict.keys())[i] != list(taut_pair_dict.keys())[j]:
                    #f.write('An item in {0} was found in {1}\n'.format(list(taut_pair_dict.keys())[i], list(taut_pair_dict.keys())[j]))
    f.close()
    # The following part has complexity O(n^2)
    # TODO: try to see if you can turn this into O(n logn)
    # This only removes the tautomers (all except the most stable one in the graph)
    for taut_list, smiles_list in taut_pair_dict.items():
        sorted_smiles = []
        for item in taut_list:
            if item in smiles_list:
                sorted_smiles.append(item)
        for index, smiles in reversed(list(enumerate(sorted_smiles))):
            # keep the highest scoring item (which is the first item, i.e. index 0 in the list)
            if index != 0:
                # Caution: the following line is very ugly
                canonical_dgvertex = dg_vertices_dict[inv_smiles_graph[rdkit_mod_smiles[sorted_smiles[0]]]]
                to_remove.update({smiles: canonical_dgvertex})

    for graph, mod_smiles in subset_mod_smiles.items():
        if mod_rdkit_smiles[mod_smiles] in to_remove:
            inEdges = [e for e in dg_vertices_dict[graph].inEdges]
            # subset only contains the graphs of the recent generation, so the list of outEdges is empty 
            for e in inEdges:
                for source in e.sources:
                    d = Derivations()
                    d.left = [source.graph]
                    d.rules = e.rules
                    d.right = [to_remove[mod_rdkit_smiles[mod_smiles]].graph]
                    b.addDerivation(d)
                    print('Added fake edge', d, 'with rule', [rule.id for rule in e.rules])
            subset.remove(graph)
            universe.remove(graph)
            print(f"Removed {graph} from subset and universe")
    return subset, universe'''

# Number of generations we want to perform
generations = 3

postSection('Final Network')
dg = DG(graphDatabase=inputGraphs)
with dg.build() as b:
    res = b.execute(strat, verbosity=2, ignoreRuleLabelTypes=True)
    subset, universe = clean_taut(dg, res)
    for gen in range(generations-1):
        res = b.execute(addSubset(subset) >> addUniverse(universe) >> strat,
                            verbosity=2, ignoreRuleLabelTypes=True)
        print('Original subset size:', len(res.subset))
        subset, universe = clean_taut(dg, res)
        print('Subset size after removal:', len(subset))
        # This step replaces the previous subset (containing tautomers) with the cleaned subset
        res = b.execute(addSubset(subset) >> addUniverse(universe))
    print('Completed')
#dg.print()
postSection('Individual Vertices')
p = GraphPrinter()
p.simpleCarbons = True
p.withColour = True
p.collapseHydrogens = True
'''for fb in forbidden:
    fb.print()
for v in dg.vertices:
    for item in forbidden:
        if v.graph.monomorphism(item) == 1:
            print(f'Found substructure {item.name} in {v.graph.name}')'''
'''for v in dg.vertices:
    v.graph.print(p)
postSection('Individual Edges')
for e in dg.edges:
    e.print()'''