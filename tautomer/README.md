# tautomer

Work on tautomer handling, kept separate from the generation pipeline.

`molvs/` tries [MolVS](https://molvs.readthedocs.io/) for tautomer canonicalisation, so that
species differing only by a tautomeric shift are not counted twice. `molvs/test.ipynb` is
the comparison that was run, next to the note that came with it.
