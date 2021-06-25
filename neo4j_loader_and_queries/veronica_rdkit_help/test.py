from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

mol = Chem.MolFromSmiles('OCC1OC(O)C(O)C(O)C1O') # glucose: https://www.researchgate.net/figure/Different-types-of-string-representations-of-the-structure-of-b-D-glucose-6_fig2_233746055

# ctrl+f for alcohol/ether: https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
hydroxyl = Chem.MolFromSmarts('[#6][OX2H]')
ether = Chem.MolFromSmarts('[OD2]([#6])[#6]')
hydroxyl_matches = mol.GetSubstructMatches(hydroxyl)
ether_matches = mol.GetSubstructMatches(ether)

print("\nHydroxyl matches:")
print(hydroxyl_matches)

print("\nEther matches:")
print(ether_matches)

# only count the entries in hydroxyl_matches for which
# the same entry does not exist in ether_matches
filtered_hydroxyl_matches = []
for hm in hydroxyl_matches:
    if not hm in ether_matches:
        filtered_hydroxyl_matches.append(hm)

print("\nHydroxyl matches without any possible ether matches:")
print(filtered_hydroxyl_matches)


from notebook_code import *


