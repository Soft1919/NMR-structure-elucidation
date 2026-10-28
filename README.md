# nmr-structure-elucidation

A Claude Code skill for solving small-molecule structures from 1D/2D NMR data, with tools that
read spectra exported as vector graphics (EMF) from Word/PowerPoint files and test candidate
structures exhaustively.

## What it does

- **Exact data extraction.** `scripts/emf.py` pulls EMF vector spectra out of .docx/.pptx;
  `scripts/spectra.py` calibrates axes automatically and gives 1D peak lists, multiplet
  splittings, integrals, 2D cross-peak lists, projections and zoomed re-plots.
- **Constitution testing.** `scripts/topology.py` tries every assignment of the observed
  signals to each candidate structure and reports the minimum number of violated HMBC (>3J)
  and COSY correlations.
- **Relative configuration.** `scripts/stereo.py` enumerates stereoisomers, removes those that
  cannot be built in 3D, averages MMFF conformers and scores Haasnoot-Altona Karplus couplings
  and NOE distances against the data.
- **Workflow.** `SKILL.md` is the procedure Claude follows; `references/checklist.md` is the
  short checklist plus config formats.

## Install

```bash
git clone https://github.com/Soft1919/NMR-structure-elucidation.git ~/.claude/skills/nmr-structure-elucidation
pip install -r ~/.claude/skills/nmr-structure-elucidation/requirements.txt
```

Then ask Claude Code to solve or check a structure from your spectra.

## Tests

```bash
python -I tests/test_karplus.py        # HLA Karplus vs. pyranose couplings (RMSD < 1 Hz)
python -I tests/test_emf.py            # EMF parser on a synthetic metafile
python -I scripts/stereo.py examples/stereo_glucose.json --nconf 20
```

## Limits

- EMF extraction assumes black axes with numeric tick labels (true for MestReNova, TopSpin and
  most plotting exports). Check the calibration residual it prints.
- Karplus and force-field results carry roughly 1-2 Hz error; use them to compare isomers.
- NMR gives relative configuration only.

## Privacy

The tools never copy embedded OLE objects (e.g. MestReNova documents), which may contain user
names, local paths and sample IDs. Do not commit source documents, spectra or work folders;
`.gitignore` excludes `work/` and common data formats.

## License

MIT
