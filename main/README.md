### Main Directory
* Each folder in this directory was meant to be for running a separate reaction (e.g. "glucose" is for glucose degradation).
* [main.py](main.py) contains methods and settings that will be used by reactions running in all subdirectories. It also contains some utility methods useful for seeing what's going on and debugging.
* [mod_to_neo4j_exporter.py](mod_to_neo4j_exporter.py) contains method(s) to create exports of reaction information that was later loaded into a [Neo4j](https://neo4j.com/) module for running graph queries (searching for autocatalytic cycles, etc.). The format of the *rels_i.txt* (where *i* is a number) files for reactions is
```
REACTION_ID MOLECULE_SMILES CONSUMPTION/CREATION    RULE_NAME
```
Here, ```REACTION_ID``` looks something like 16_0, where 16 is the id of the hyperedge that connects the molecules (in MOD) and the subscript 0 is an index (running from 0 to n-1) to specify a given rule associated with the hyperedge. For instance, a particular set of edges is associated with amidine hydrolysis as well as some sort of tautomerisation reaction, they could each be assigned a separate index in the reaction id.
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

##TODO:
* explain the format of glucose_degradation_output.txt
* the (redundant, as of now) tautomer cleaning methods.