# Reaction Space Explorer

## TODO:
* ~~Get rid of molecules containing unstable substructures (three or four membered rings)~~
    * Currently only a few substructures are forbidden, make the library of unstable/impossible substructures more complete.
* Fix the reaction rules that aren't currently invertible for some reason.
    * Check the ones with wildcards.
* Write more reaction rules (to account for missing products).
* ~~Compare with mass spectra at the end of each generation and see what % of masses are matching.~~
    * ~~Perhaps represent this info graphically.~~
* ~~See which step in the tautomer cleaning process is the most inefficient.~~
    * The slowest part is enumerating all possible tautomers for the given molecules. Takes more time to enumerate the possibilties than it does to produce the molecules themselves.
