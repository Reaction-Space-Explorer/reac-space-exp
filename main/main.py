import os
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from molvs import tautomer

include(os.path.abspath(os.path.join('..', 'rules/all.py')))
postChapter('Alkaline Glucose Degradation')

# starting molecule
glucose = smiles('OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O', name='Glucose')

# Forbidden substructures
# Three and four membeered rings are unstable
forbidden = [smiles('C1CC1', name='cyclopropane'), smiles('C1CCC1', name = 'cyclobutane')]

# make sure these don't get passed as an input
for fb in forbidden:
    inputGraphs.remove(fb)

strat = (
    addSubset(inputGraphs)
    >> rightPredicate[
        lambda derivation: all (g.exactMass <= 500
         and (g.monomorphism(fb) == 0 for fb in forbidden) for g in derivation.right)
    ] (repeat[1](inputRules))
)

# molVs smiles might not be identical to rdkit so let's use RDKit instead of molVs
enum = rdMolStandardize.TautomerEnumerator()
#enum = tautomer.TautomerEnumerator()

# What if none of the existing tautomers is the most stable one?
# We keep the most stable out of the ones in the network since we can't modify graphs in MOD
# enumerate all possible tautomers for each molecule
#   find if tautomeric pairs exist for each of them
#       if true: sort the tautomers and remove the ones with the lowest score
# items are to be removed from both subset and universe
# return the cleaned tautomers once done.

# The code became really ugly while optimizing for performance using dictionaries
def clean_taut(dg, dg_execute):
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
    return subset, universe

# Number of generations we want to perform
generations = 2

dg = DG(graphDatabase=inputGraphs)
with dg.build() as b:
    res = b.execute(strat, verbosity=2, ignoreRuleLabelTypes=True)
    subset, universe = clean_taut(dg, res)
    for gen in range(generations-1):
        print(f'Graph has total {dg.numVertices} vertices')
        res = b.execute(addSubset(subset) >> addUniverse(universe) >> strat,
                            verbosity=2, ignoreRuleLabelTypes=True)
        print('Original subset size:', len(res.subset))
        subset, universe = clean_taut(dg, res)
        print('Size after removal:', len(subset))
        # This step replaces the previous subset (containing tautomers) with the cleaned subset
        res = b.execute(addSubset(subset) >> addUniverse(universe))
    print('Completed')
#dg.print()