# Reaction Space Explorer
An open source cheminformatics workflow to simulate chemical reaction networks important in prebiotic chemistry and to discover autocatalytic loops. Available open source under the BSD-3 Clause license.

The platform uses:
* [MØD](https://github.com/jakobandersen/mod) version [v.0.11.0](https://github.com/jakobandersen/mod/releases/tag/v0.11.0).2 for network generation.
* [RDKit](https://anaconda.org/rdkit/rdkit) (not compatible with Python 3.8.x as of now, we used a separate conda environment with Python 3.7.9).
* [Neo4j](https://neo4j.com/) for graph queries.

Figures were plotted using [matplotlib](https://matplotlib.org/).
### Note:
The FT-ICR-MS data has some rows with inconsistent number of columns. For the purpose of plotting the mass spectra, we extracted columns of our interest.

## Publications
To-be-updated.
## TODO:

* ~Try to match all the possible structures in the Y&M paper.~
    * 96% structures have matched for glucose degradation, 100% for formose.
* For the future take into account Kinetics for the molecules in the network (https://rmg.mit.edu/) 
