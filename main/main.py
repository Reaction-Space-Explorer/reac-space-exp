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
    for smiles in to_remove:
        # The DGVertex associated with the molecule to be removed
        dg_vertex = mod_dgverts[graphs[mod_subset_smiles.index(smiles)]]
        for e in dg_vertex.inEdges:
            for source in e.sources:
                    d = Derivations()
                    d.left = [source.graph]
                    d.rules = e.rules
                    d.right = [graphs[mod_subset_smiles.index(smiles)]]
                    b.addDerivation(d)
        subset.remove(graphs[mod_subset_smiles.index(smiles)])
        universe.remove(graphs[mod_subset_smiles.index(smiles)])
    return subset, universe

# Number of generations we want to perform
generations = 2

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