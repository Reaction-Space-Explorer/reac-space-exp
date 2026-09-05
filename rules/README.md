# rules

The reaction rule library: 54 files, GML rules for MØD, one transformation or a small family
per file, named for the chemistry — `aldolCondensation.py`, `ketoEnolisation.py`,
`retroAldol.py`.

- `all.py` — calls every rule file, so loading it loads the library
- `common.py` — helpers the rule files share

To add a rule, write it in its own file following an existing one and add the call to
`all.py`. A visual list is at [`data/All_Rules.pdf`](../data/All_Rules.pdf).
