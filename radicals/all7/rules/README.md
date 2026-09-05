# rules

The 19 rules this run loads, one transformation per file, kept separate from the driver so
`all7.py` stays readable. They are GML rules as in the closed-shell library at
[`rules/`](../../../rules), but acting on radicals: the chemistry here is oxygen, ozone and
the HOx radicals, and each rule's `ruleID` names the step it encodes.

| | Rule | |
|---|---|---|
| `foto-1.py`, `foto-2.py` | `foto O2` | O₂ photolysed to two oxygen atoms |
| `foto-3.py` | `foto O3` | O₃ photolysed |
| `abstforox.py`, `ozone-rad1.py` | `O3 attack by O` | Ozone attacked by an oxygen atom |
| `ozone-rad2.py` | `O3 attack by .O-X, X != O` | Ozone attacked by an oxygen-centred radical |
| `ozone-rad3.py` | `O3 attack by .H` | Ozone attacked by a hydrogen atom |
| `oxygen-rad1.py`, `oxygen-rad2.py` | `O attack by .O-H`, `.O-OH` | Atomic oxygen attacked by hydroxyl and by hydroperoxyl |
| `h-abstr-rad1.py` | `.O2H attack by H.` | Hydroperoxyl attacked by a hydrogen atom |
| `h-abstr-rad2.py`, `h-abstr-rad3.py` | `.O2H`, `.OH attack by .O-X` | Hydroperoxyl and hydroxyl attacked by an oxygen-centred radical |
| `termination-rad1.py`, `termo-hy-rad*.py`, `termo-oxygen-rad*.py` | `termination`, `termo …` | Two radicals combining, ending the chain |

`foto` is photo and `termo` termination; the names are Spanish, as elsewhere in
[`radicals/`](../..). `abstforox.py` and `ozone-rad1.py` carry the same rule, and
`termo-hy-rad1.py` and `termo-oxygen-rad4.py` likewise.
