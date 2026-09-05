# radicals

Radical chemistry networks, generated with MØD like the rest but kept separate: the rules
here act on radical species rather than closed-shell molecules.

```bash
mod -f all6.py
```

MØD appends to its output rather than replacing it, so a previous run has to be cleared
first, or the new network will be generated on top of the old one:

```bash
rm -r Rule* Density* dgdumfile-* Masses* Stage*
mod clean
```

`Etapa0.txt` then has to be recreated with the starting molecules.

| | |
|---|---|
| `all6.py` | The driver. Includes the rule files and plots how molecular mass evolves over generations |
| `common.py` | Helpers the rule files share, included by each |
| `Etapa0.txt` | Starting molecules, read at generation 0 |
| `OHRadAlcanes.py` | Hydroxyl radical abstracting from an alkane |
| `OHRadHM.py`, `OHRadO=M.py` | Hydroxyl radical against the other substrate classes |
| `doublebond-attack-Mrad.py` | Radical addition across a double bond |
| `sustract-proton.py` | Proton abstraction |
| `termination.py` | Radical recombination, ending a chain |
| `all7/` | A later run, see its own README |

Several files and comments are in Spanish; `Etapa` is stage, `sustract` is abstract.
