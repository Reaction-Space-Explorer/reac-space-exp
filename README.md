# Reaction Space Explorer

## TODO:
* Try to match all the possible structures in the Y&M paper.
* ~~For some reason, comparison of molecules in the SDF that contain Y&M's reported structures isn't working correctly, fix that.~~
* ~~Get rid of molecules containing unstable substructures (three or four membered rings)~~
    * ~~Currently only a few substructures are forbidden, make the library of unstable/impossible substructures more complete~~. Made all three and four membered rings, atoms with two double bonds forbidden using wildcards.
* ~~Fix the reaction rules that aren't currently invertible for some reason~~.
    * ~~Check the ones with wildcards~~. Rules with ```constrainAdj[]``` will not be directly invertible.
* Write more reaction rules (to account for missing products).
    * ~~Hemiacetal formation for 6 membered rings~~.
    * ~~Aldol/Retro-aldol~~
    * ~~Elimination~~/Hydration (one hydration rule left).
    * Ester Hydrolysis (inverse)
    * ~~Keto-enol migration~~
        * Maybe avoid enols by making it go two steps (requires adjacent -OH available).
* ~~Compare with mass spectra at the end of each generation and see what % of masses are matching.~~
    * ~~Perhaps represent this info graphically.~~
* ~~See which step in the tautomer cleaning process is the most inefficient.~~
    * The slowest part is enumerating all possible tautomers for the given molecules. Takes more time to enumerate the possibilties than it does to produce the molecules themselves.
