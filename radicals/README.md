# radicals

Radical networks, generated with MØD like the rest but with rules acting on radicals rather
than closed-shell molecules. Run with `mod -f all6.py`.

MØD appends to its output, so clear a previous run first, then recreate `Etapa0.txt` with the
starting molecules:

```bash
rm -r Rule* Density* dgdumfile-* Masses* Stage*
mod clean
```

- `all6.py` — the driver; includes the rule files and plots how mass evolves over generations
- `common.py` — helpers included by each rule file
- `Etapa0.txt` — starting molecules
- `OHRadAlcanes.py`, `OHRadHM.py`, `OHRadO=M.py` — hydroxyl radical against each substrate class
- `doublebond-attack-Mrad.py` — radical addition across a double bond
- `sustract-proton.py` — proton abstraction
- `termination.py` — radical recombination
- `all7/` — a later run

Names are Spanish: `Etapa` is stage, `sustract` abstract.
