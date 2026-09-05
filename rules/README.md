# rules

The reaction rule library, 54 files, written as GML graph transformation rules for MØD.

Each `.py` file holds one transformation, or a small family of them, named for the chemistry
it encodes: `aldolCondensation.py`, `ketoEnolisation.py`, `retroAldol.py` and so on.

- `all.py` calls each rule file, so loading it loads the whole library.
- `common.py` holds the helpers the rule files share, and `all.py` imports them.

To add a rule, write it in its own file following an existing one, and add the call to
`all.py`. The MØD [documentation](https://jakobandersen.github.io/mod/) has worked examples
of rule syntax. A visual list of the library is at [`data/All_Rules.pdf`](../data/All_Rules.pdf).
