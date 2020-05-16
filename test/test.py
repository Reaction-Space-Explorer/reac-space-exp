# MOD Example code
# Normal printing to the terminal:
print("Hello world")
# Make some headers in the summary:
postChapter("Hello")
postSection("World")
# Load a moleucle from a SMILES string:
mol = smiles("Cn1cnc2c1c(=O)n(c(=O)n2C)C", name="Caffeine")
# Put a visualisation of the molecule in the summary:
mol.print()
