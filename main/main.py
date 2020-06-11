import os
import rdkit.Chem as chem

from molvs import tautomer
from molvs import standardize
from molvs import standardize_smiles

include(os.path.abspath(os.path.join('..', 'rules/all.py')))
postChapter('Alkaline Glucose Degradation')

# Create a dictionary of smiles strings (automatically avoids duplicates)
# with the value being the generation in which the molecule was produced
molecule_smiles = {}
cleaned_smiles = {}

# starting molecule
glucose = smiles('OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O', name='Glucose')

print(type(inputGraphs))
# assign zeroth generation to input molecules
molecule_smiles.update({glucose.smiles : 0})

# TODO: This only removes the tautomers from a list of smiles sent externally from MOD
# We need to cause a "side effect" on the derivation graph to safely remove the useless nodes
# The properties of the DG could be dumped and reloaded to achieve this.
generation = 1

strat = (
    addSubset(inputGraphs)
    >> rightPredicate[
        lambda derivation: all (g.exactMass <= 500 for g in derivation.right)
    ] (repeat[2](inputRules))
)
dg = DG(graphDatabase=inputGraphs)
dg.build().execute(strat, ignoreRuleLabelTypes=True)

for v in dg.vertices:
    molecule_smiles.update({v.graph.smiles: generation})

generation += 1

# Clean tautomers
canon = tautomer.TautomerCanonicalizer()
for smiles, gen in molecule_smiles.items():
    canonical = canon.canonicalize(chem.MolFromSmiles(smiles))
    canonical_taut = chem.MolToSmiles(canonical)
    # if the smiles is already in the dictionary, it won't be updated
    # TODO: think of what would happen if the same molecule popped in two different generations
    cleaned_smiles.update({canonical_taut: gen})

print(f'Number of tautomers removed: {len(molecule_smiles.keys()) - len(cleaned_smiles.keys())}')

