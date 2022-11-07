# Reaction Space Explorer
An open source cheminformatics workflow to simulate chemical reaction networks important in prebiotic chemistry and to discover autocatalytic loops. Available open source under the BSD-3 Clause license.

The platform uses:
* [MØD](https://github.com/jakobandersen/mod)
    * Note: Most of the reaction networks were generated with MØD versions v.0.9.0 to [v.0.11.0](https://github.com/jakobandersen/mod/releases/tag/v0.11.0). However, the code (and any `.dg` output files that may have used older formats) are still supported. We have provided `.dg` files with newer formats that come with reaction rules self contained in them. Note that for reproducing reactions from scratch you will still need to load rules from the [library we compiled](rules/). The user is free to add or remove any rules as per their choice.
    * Although MØD installation has been described on its documentation pages, some people have benefitted from a [small tutorial](Setting_up_MOD.md) we ourselves wrote.
* [RDKit](https://anaconda.org/rdkit/rdkit) (not compatible with Python 3.8.x as of now, we used a separate conda environment with Python 3.7.9).
* [Neo4j](https://neo4j.com/) for graph queries.
* [Gephi](https://gephi.org/) for network visualization.

Figures were plotted using [matplotlib](https://matplotlib.org/), [seaborn](https://seaborn.pydata.org/), and [skunk](https://github.com/whitead/skunk).

## Organization of this repository
### Scripts to run reactions
The [main](main/) folder contains subfolders for each reaction (e.g. [glucose degradation](main/glucose/), and [formose](main/formose/)) that we have studied. The output produced by running this pipeline has been placed in subfolders inside those.

### Available output files
For reactions with published or submitted manuscripts, we have provided.
* A tab-separated table containing SMILES of each species and the generation in which it was produced.
* The `.dg` file that was dumped via MØD itself, which can be used to load the full network into MØD without the need of generating it again using the [DG.load](https://jakobandersen.github.io/mod/pymod/dg/DG.html?highlight=dg%20load#mod.DG.load)() method built into it.
* A custom format (tab-separated) output that can be loaded in Neo4J for network queries, or Gephi for the purpose of visualization.
* A table of rules applied by generation.

The methods for generating these are contained in the code scripts

See the README description in [main](main/) for an overview of the formats.

### Library of reaction rules
All rules are placed in the [rules](rules/) folder. They are written in GML format, mostly compartmentalized in separate `.py` files. For convenience purposes, we used some tricks to make things modular. The `all.py` file calls each `.py` file in the folder. Additional rules that you'd like to load (*if* using our methods) would The same file also imports utility methods from `common.py`. Examples that demonstrate how to write rules can be found on the [documentation pages of MØD](https://jakobandersen.github.io/mod/).

### Complementary data
The [data](data/) folder contains a list of forbidden substructures that we utilize for post-generation filtering (discussed in Arya *et al.* 2022), a visual list of all reaction rules (which may be slightly dated), table of thermochemical data, among other things.

### Autocatalysis Search and Figures
The [imperative_loader_and_queries](imperative_loader_and_queries/) folder contains the files, and instructions for spontaneous autocatalytic loop search, among other things.


### Running a file using MØD
The way we ran the files were using the terminal 
```bash
mod -f glucose_degradation.py
```
You will notice there are some weird calls to methods that you may think have not even been imported, such as ```include("main.py")```. When a python file is run via MØD, it auto imports packages and methods in *libPyMØD*.

## Publications
* Arya, A., Ray, J., Sharma, S., Cruz Simbron, R., Lozano, A., Smith, H. B., ...Cleaves, H. J. (2022). An open source computational workflow for the discovery of autocatalytic networks in abiotic reactions. _Chemical Science_ 13, 4838–4853. https://doi.org/10.1039/d2sc00256f
    * For this, the reaction studied was [glucose degradation](main/glucose), and the [output files](main/glucose/output) for it can be found in the relevant folder.
* Sharma, S.; Arya, A.; Cruz, R.; Cleaves II, H.J. Automated Exploration of Prebiotic Chemical Reaction Space: Progress and Perspectives. _Life_ 2021, 11, 1140. https://doi.org/10.3390/life11111140

Stay tuned for updates..

