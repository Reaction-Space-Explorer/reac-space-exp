### Main Directory
* Each folder in this directory was meant to be for running a separate reaction (e.g. "glucose" is for glucose degradation).
* [main.py](main.py) contains methods and settings that will be used by reactions running in all subdirectories. It also contains some utility methods useful for seeing what's going on and debugging.
* [mod_to_neo4j_exporter.py](mod_to_neo4j_exporter.py) contains method(s) to export networks from the MØD pipeline that we  later loaded into a [Neo4j](https://neo4j.com/) MØDule for running graph queries (searching for autocatalytic cycles, etc.).

### Custom dump file formats
**Species by generation**:

It was important to keep track of which species was created in which generation. ```.txt``` files with names such as ```reaction_name_output.txt``` contain exactly this. For instance, here's an excerpt from ```glucose_degradation_output.txt```

```
G1	C(C(C1C(C(C(O)O1)O)O)O)O
G1	C(C1C(C(C(C(O)O1)O)O)O)O
G1	C1C(C(C(C(C(O)O1)O)O)O)O
G2	C1(C(C(C(C1CO)O)O)O)=O
G2	C1(C(C(C(C(C1)O)O)O)O)=O
G2	C(C1C(C(C(C1O)O)O)O)=O
```
Here, G1 means the corresponding SMILES (unique for each species, canonicalized by MØD) first appeared in generation 1. ```G0``` would include just the initial reactants.

**Neo4j Imports**:

The neo4j exporter creates two kinds of files, one which shows *relations* between species (reactions) and others just show which species exist. The format of the *rels_i.txt* (where *i* represents the *generation*) files for reactions is
```
REACTION_ID MOLECULE_SMILES CONSUMPTION/CREATION    RULE_NAME
```
A description of the format is as follows:

A given set of products connected to a set of edges could be formed not by just one reaction rule but multiple rules (which involve the same reactants/products). To account for that, we added a subscript to the reaction IDs where each subscript value corresponds to a different rule. Say, a given reaction edge has 'n' rules associated with it. We create separate entries for the reaction with subscripts running from 0 to n-1 

Here, ```REACTION_ID``` looks something like 16_0, where 16 is the id of the hyperedge that connects the molecules (in MØD) and the subscript 0 means it's the first rule associated with that hyperedge. It was important to make a distinction because particular set of edges could associated with amidine hydrolysis as well as some sort of tautomerisation reaction, they could each be assigned a separate subscript index for the same reaction id.
```MOLECULE_SMILES``` is simply the SMILES string associated with that molecule.
```CONSUMPTION/CREATION``` is -1 if the molecule is a *reactant* and is 1 if it's a *product*.
```RULE_NAME``` is simply the name of the reaction rule which produced it.

An example is
```
6_0	C#N	-1	HCN Addition to Nitriles
6_0	C#N	-1	HCN Addition to Nitriles
6_0	C(C=N)#N	1	HCN Addition to Nitriles
```
which represents ```C#N + C#N -> C(C=N)#N``` happening.

**Note**: As of now, the contents of ```rels_n.txt``` are cumulative in the sense that a file, say ```rels_3.txt``` would contain the content of ```rels_2.txt``` and ```rels_1.txt``` as well. So ```rels_3.txt``` includes reactions of G1 + G2 + G3 combined. Therefore, care should be taken when trying to obtain just the reactions of G3.

**Precaution**: every subsequent run, these ```.txt``` dump files will not be deleted and rewritten on their own, instead additional lines will be appended at their end. So one needs to delete old files manually.

## Tautomer Cleanup
We experimented with different techniques (using RDKit/MolVS as well as Ambit-tautomer) to remove redundant tautomers but none of them were fully satisfactory. Some of the (redundant) code which we no longer employ in our workflow is still here (e.g. [clean_tautomers.py](clean_tautomers.py), [ambit_tautomer.py](ambit_tautomer.py)). Ambit is a Java library (dependent on CDK) which, to our workflow, could be pipelined using ```jpype```. See [Using_Ambit.md](Using_Ambit.md) for notes that I wrote in mid-2020 if they seem useful to you.
## TODO:
* see if there's anything else useful to explain

## Future Code Improvements
* The way I wrote the "cheap hax" to toggle the `Cannizarro 2` rule by setting `with_formaldehyde=True` or `=False` was terrible. This boolean was defined initially in `../rules/cannizarro2.py`, to decide whether to use glucose as a "fixed aldehyde" or to use formaldehyde instead. Things can easily go wrong because of how this boolean has been accessed/mutated in other scripts like `main.py` or the individual reaction scripts like `glucose_degradation.py`.
* The way these scripts have been run up to now is by performing `mod -f filename.py` --- that still provides most of Python functionality but I have noticed conflicts with MOD's way of dealing with package directories and Python's native way. In newer versions, importing `*` from `libpymod` directly could reduce that.
    * A more proper way of dealing with packages would then be possible.
* A lot of what I put in `main.py` (back in 2020) are 'utility' methods. These could go in their own separate file.
* Several lengthy methods could be modularised better for the sake of readability. Also, I learnt much cleaner ways to do certain tasks later on that I haven't updated everywhere. (These don't affect the function of the code itself, but make it more readable.)