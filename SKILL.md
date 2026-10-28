---
name: nmr-structure-elucidation
description: Solve small-molecule structure elucidation problems from 1D/2D NMR (1H, 13C, DEPT, HSQC, COSY, HMBC, NOESY) plus MS/IR, including spectra embedded as vector images in Word/PowerPoint files. Use when the user gives spectra and a formula and asks for the structure, the relative configuration, a signal assignment, or a check of a proposed/published structure.
---

# NMR structure elucidation

Work from the spectra. Do not look up the answer or try to identify the compound in the
literature unless the user asks; a published structure can be wrong, and the deliverable is a
structure that the data support. State clearly what is proven, what is inferred, and what is open.

## 0. Get exact numbers

- Spectra inside .docx/.pptx are often EMF vector graphics. Extract them instead of reading
  pictures: `python scripts/emf.py INPUT.docx --out work/emf`, then
  `python scripts/spectra.py peaks1d|mult|integ|peaks2d|proj|zoom ...` (see `--help`).
- Check the printed calibration residuals. Note the digital resolution of each trace
  (decimated traces merge close peaks and under-integrate sharp singlets).
- 2D plots usually carry their own 1D projections (`proj`): a second, independent peak list
  and a check that 1D and 2D data come from the same sample.
- Raster-only spectra: read them, but quote the reading uncertainty.

## 1. Remove what is not the compound

Before interpreting, list and exclude: solvent and labelled impurities; residual CH2Cl2 / grease /
water; HMBC one-bond artefacts (pairs symmetric about a 1H shift, split by 1J(CH) ~125-180 Hz,
sitting on that proton's own HSQC carbon); DEPT peaks with no 13C counterpart; NOESY cross-peaks
with the same sign as the diagonal (exchange/COSY artefacts in small molecules).

## 2. Inventory

1. Formula -> degrees of unsaturation; nitrogen rule; IR functional groups.
2. Every carbon: 13C + DEPT-135 sign + HSQC protons -> C / CH / CH2 / CH3. Overlapped carbons
   show up as mixed DEPT signs or extra HSQC protons; the carbon count must equal the formula.
3. Count C-bound H; the difference to the formula = OH/NH. Distribute heteroatoms
   (shift ranges) and the remaining unsaturation (C=X vs rings).

## 3. Constitution

1. COSY + HSQC -> spin fragments; note where each fragment ends (quaternary C, heteroatom).
2. HMBC joins fragments. Use only 2J/3J; methyl singlets are the most reliable anchors.
3. Use valence bookkeeping: when only two atoms have open valences, they must be bonded.
4. Write down every competing constitution that fits the formula and class of compound.
   Test all of them with `scripts/topology.py` (exhaustive assignment search; reports the
   minimum number of violated HMBC/COSY correlations). Zero violations is required for the
   proposed structure; every violation of a rival must be a real, checked correlation.
5. Look for correlations that must exist for a candidate and are absent (e.g. an expected
   large vicinal coupling or COSY cross-peak). Absences are often decisive.

## 4. Relative configuration

1. Qualitative first: large 3J (>8 Hz) = near-antiperiplanar; NOE contacts within ~3 A.
   Rigid polycycles: decide which face each substituent occupies.
2. Then exhaustive: `scripts/stereo.py CONFIG.json` enumerates stereoisomers, discards
   non-embeddable ones, Boltzmann-averages MMFF conformers and scores HLA-Karplus J values and
   r^-6 distances against the data. Check robustness: re-weight NOE vs J, drop one observable
   at a time; report which observable decides each centre.
3. NMR gives relative configuration only. Never state an absolute configuration from NMR.

## 5. Report

- Final structure (SMILES with stereo + a drawing), full assignment table
  (dC, type, dH (mult, J), key HMBC, key NOE), and the argument chain step by step.
- A short "not proven" list: interchangeable assignments, overlapped signals, configurations
  that rest on a single observable, absolute configuration.
- If a reference/published structure disagrees, list the specific correlations it violates
  and the experiment that would settle it. Do not adopt it silently.

## Pitfalls

- Drawing: build 2D coordinates yourself and re-read the molblock to confirm stereo; removing
  H atoms one by one can corrupt parity (use RemoveHs). Label substituent carbons too.
- Karplus values carry ~1-2 Hz model error; compare isomers with one consistent method.
- Do not trust a single weak cross-peak, and do not ignore a strong one.
- Never redistribute the source documents or embedded OLE/MestReNova objects: they can
  contain user names, file paths and sample IDs. Strip them before sharing anything.

See `references/checklist.md` for a printable checklist and the config formats.
