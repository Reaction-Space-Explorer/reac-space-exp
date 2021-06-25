# Reaction Space Explorer
An open source cheminformatics workflow to simulate chemical reaction networks important in prebiotic chemistry and to discover autocatalytic loops. Available open source under the BSD-3 Clause license.

The platform uses:
* [MØD](https://github.com/jakobandersen/mod) version [v.0.11.0](https://github.com/jakobandersen/mod/releases/tag/v0.11.0).2 for network generation. (More recent versions are also supported.)
* [RDKit](https://anaconda.org/rdkit/rdkit) (not compatible with Python 3.8.x as of now, we used a separate conda environment with Python 3.7.9).
* [Neo4j](https://neo4j.com/) for graph queries.

Figures were plotted using [matplotlib](https://matplotlib.org/).

### Running a file using MØD
The way we ran the files were using the terminal 
```bash
mod -f glucose_degradation.py
```
You will notice there are some weird calls to methods that you may think have not even been imported, such as ```include("main.py")```. When a python file is run via MØD, it auto imports packages and methods in *libPyMØD*.

## Publications
To-be-updated.
## TODO:

<<<<<<< Updated upstream
* ~Try to match all the possible structures in the Y&M paper.~
    * 96% structures have matched for glucose degradation, 100% for formose.
* For the future take into account Kinetics for the molecules in the network (https://rmg.mit.edu/) 
=======
>>>>>>> Stashed changes
