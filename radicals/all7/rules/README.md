# rules

The 19 rules this run loads, one per file, GML as in [`rules/`](../../../rules) but acting on
radicals. The chemistry is oxygen, ozone and the HOx radicals; each file's `ruleID` names the
step.

- `foto-1.py`, `foto-2.py` — `foto O2`, O₂ photolysed to two oxygen atoms
- `foto-3.py` — `foto O3`
- `abstforox.py`, `ozone-rad1.py` — `O3 attack by O`
- `ozone-rad2.py`, `ozone-rad3.py` — ozone attacked by an oxygen-centred radical, and by H
- `oxygen-rad1.py`, `oxygen-rad2.py` — atomic oxygen attacked by `.OH` and by `.O-OH`
- `h-abstr-rad1.py` — `.O2H attack by H.`
- `h-abstr-rad2.py`, `h-abstr-rad3.py` — `.O2H` and `.OH` attacked by an oxygen-centred radical
- `termination-rad1.py`, `termo-hy-rad*.py`, `termo-oxygen-rad*.py` — two radicals combining

`abstforox.py` and `ozone-rad1.py` hold the same rule, as do `termo-hy-rad1.py` and
`termo-oxygen-rad4.py`. `foto` is photo, `termo` termination.
