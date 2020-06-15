import os
import mod
from mod import *
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from molvs import tautomer

include(os.path.abspath(os.path.join('..', 'rules/all.py')))
postChapter('Alkaline Glucose Degradation')

# Create a dictionary of smiles strings (automatically avoids duplicates)
# with the value being the generation in which the molecule was produced
molecule_smiles = {}
cleaned_smiles = {}

# starting molecule
glucose = smiles('OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O', name='Glucose')

# assign zeroth generation to input molecules
molecule_smiles.update({glucose.smiles : 0})

# TODO: This only removes the tautomers from a list of smiles sent externally from MOD
# We need to cause a "side effect" on the derivation graph to safely remove the useless nodes
# The properties of the DG could be dumped and reloaded to achieve this.

strat = (
    addSubset(inputGraphs)
    >> rightPredicate[
        lambda derivation: all (g.exactMass <= 500 for g in derivation.right)
    ] (repeat[1](inputRules))
)

canon = tautomer.TautomerCanonicalizer()
# Let's use RDKit instead of molVs for now
enum = rdMolStandardize.TautomerEnumerator()
# What if none of the existing tautomers is the most stable one? Do keep one based on a score?

# enumerate all possible tautomers for each molecule
#   find if tautomeric pairs exist for each of them
#       if true: sort the tautomers and remove the ones with the lowest score
# items are to be removed from both subset and universe
# return the cleaned tautomers once done.
def clean_taut(dg_execute):
    subset = dg_execute.subset
    universe = dg_execute.universe
    # A dictionary of Graph objects and the corresponding SMILES string produced by MOD
    subset_mod_smiles = {g: g.smiles for g in subset}
    # dictionary value will later be updated to its canonical smiles
    mod_rdkit_smiles = {g.smiles: None for g in subset}
    # This dictionary will contain {canonical rdkit smiles: [list of taut smiles]}
    smiles_tauts_dict = {}
    # list of canonical smiles of the isomers that need to be removed
    to_remove = []
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
            #to_remove.append(smiles)
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
                to_remove.append(smiles)
                #print(f'Removing {smiles} from tautomer set {sorted_smiles}')
    #TODO: add fake edges to account for the edges that were washed away with the removed nodes
    # This needs to be done before the graphs are removed.
    for graph, mod_smiles in subset_mod_smiles.items():
        if mod_rdkit_smiles[mod_smiles] in to_remove:
            subset.remove(graph)
            universe.remove(graph)
            print(f"Removed {graph} from subset and universe")
    return subset, universe

# Number of generations we want to perform
generations = 3

dg = DG(graphDatabase=inputGraphs)
with dg.build() as b:
    res = b.execute(strat, verbosity=2, ignoreRuleLabelTypes=True)
#    for g in dg.vertices:
#        g.graph.print()
    subset, universe = clean_taut(res)
    for gen in range(generations-1):
        print(f'Graph has total {dg.numVertices} vertices')
        # Now clean tautomers
        smiles_list = [g.smiles for g in subset]
        res = b.execute(addSubset(subset) >> addUniverse(universe) >> strat,
                            verbosity=2, ignoreRuleLabelTypes=True)
        print('Original subset size:', len(res.subset))
        subset, universe = clean_taut(res)
        print('Size after removal:', len(subset))
        # This step replaces the previous subset (containing tautomers) with the cleaned subset
        res = b.execute(addSubset(subset) >> addUniverse(universe))
        print('Created final version of DG')
    print('Completed')
#dg.print()