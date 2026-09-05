# Stoichiometric network analysis

Stoichiometric analysis of the networks generated in this repository, by **Alejandro Lozano**
([@alejandro-lozano-dev](https://github.com/alejandro-lozano-dev)). It reads a network as its
stoichiometric matrix and works from the four fundamental subspaces of that matrix; Gilbert
Strang's [18.06 lectures](https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/) are
the clearest introduction to them.

| | |
|---|---|
| `glucose_3rd_generation/` | Glucose degradation network at the third generation |
| `glucose_5th_generation/` | The same at the fifth generation, with the species matched and missed against the reference list, and a MØD against RDKit tautomer comparison |
| `radical_network/` | The radical network, same treatment |
| `glucose_degradation_output.txt` | The network the analyses read |

Each of the three holds the same four files: `Data_cleaning.ipynb` builds the matrix from the
network output and pickles it into `reactants` and `products`, and `Ematrix.m` finds the
extreme currents from it, leaving `test.mat`. `Ematrix.m` is Guy Schmitz's 2008 program, not
written here; each folder's README says so.

Imported from the `Stoichiometric-Network-Analysis` repository with its history. Two sets of
files from that repository are not here, both being third-party course material rather than
work on these networks: the lecture PDFs from Bernhard Palsson's *Dynamic States of
Biological Networks*, and the worked homework for that course. The
[MassPy documentation](https://masspy.readthedocs.io/en/latest/notebooks/SB2_textbook/SB2-Chapter-1.html)
carries the textbook notebooks openly.
