# Checklist

## Data
- [ ] Formula, DBE, nitrogen rule, IR bands
- [ ] Spectra extracted numerically; calibration residuals noted; digital resolution noted
- [ ] 1D and 2D data shown to belong to the same sample (projection peaks agree)
- [ ] Impurities / solvent / one-bond HMBC artefacts / DEPT extras excluded and listed

## Constitution
- [ ] Every carbon typed; count matches formula; overlaps resolved
- [ ] C-bound H count vs formula -> number of OH/NH
- [ ] Spin fragments from COSY, with their ends
- [ ] HMBC links, each checked as 2J/3J
- [ ] Rival constitutions listed and tested with topology.py
- [ ] Expected-but-absent correlations checked for each rival

## Configuration
- [ ] Large J and NOE contacts interpreted qualitatively
- [ ] stereo.py run; robustness checked (weights, leave-one-out)
- [ ] Deciding observable named for each stereocentre

## Report
- [ ] Structure + assignment table + argument chain
- [ ] "Not proven" list
- [ ] No private metadata in anything shared

## topology.py config (example skeleton)

```json
{
  "candidates": {"A": "[CH:1]1...", "B": "[CH:1]1..."},
  "groups": {
    "CH":  {"signals": ["c1", "c2"], "atoms": {"A": [9, 13], "B": [9, 12]}},
    "CH2": {"signals": ["c3", "c4"], "atoms": [14, 2]}
  },
  "fixed": {"c5": 16},
  "aliases": {"h_ov": ["c6", "c7"]},
  "hmbc": [["c3", ["c5"]], ["c1", ["c6", "c7"]]],
  "cosy": [["c2", "c3"]]
}
```
Signals are free-form ids; HMBC entries are [proton-bearing carbon, [possible target carbons]].
`aliases` handles protons that overlap and could sit on either listed carbon.

## stereo.py config

See `examples/stereo_glucose.json`: `smiles` (atom-mapped, open centres unassigned), `j`
[[mapA, mapB, J]], optional `j_sets` [[mapA, [partners], [J values]]], `noe` [[mapA, mapB]],
optional `no_noe`.
