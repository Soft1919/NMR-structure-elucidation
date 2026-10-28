"""Rank stereoisomers of a fixed constitution against measured 3J(H,H) and NOE contacts.

    python stereo.py CONFIG.json [--nconf 40] [--workers 8] [--top 10]

Pipeline
  1. Enumerate stereoisomers of `smiles` (atom-mapped, may carry fixed stereo tags).
     Enantiomer pairs are collapsed: NMR cannot tell them apart.
  2. For each isomer: ETKDG embedding (isomers that cannot be embedded are geometrically
     impossible, typical for cages), MMFF94s optimisation, Boltzmann weights at 298 K.
  3. Observables: 3J by the HLA Karplus equation (karplus.py, validated on pyranoses) and
     r^-6 averaged H-H distances.
  4. Score = sum((Jcalc - Jexp)^2) + w_noe * sum(max(0, r - r_max)^2) for observed NOEs
            + w_noe * sum(max(0, r_min - r)^2) for NOEs reported as clearly absent.

CONFIG keys
  smiles   : atom-mapped SMILES (heavy-atom maps used below)
  j        : [[mapA, mapB, J_Hz], ...]    vicinal coupling between the H on A and the H on B.
             For CH2 partners the closest predicted value is used unless sorted=true.
  j_sets   : [[mapA, [mapB1, mapB2, ...], [J1, J2, ...]], ...] all couplings of one proton,
             compared after sorting (use when you do not know which partner gives which J).
  noe      : [[mapA, mapB], ...]          observed contacts (any H on each atom)
  no_noe   : [[mapA, mapB], ...]          contacts confirmed absent (use sparingly)
  r_max (3.2), r_min (3.0), w_noe (4), max_isomers (512)
"""
import argparse, json, os, sys
import numpy as np
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms as T
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
import karplus

CFG = None


def isomers(smi, limit):
    base = Chem.MolFromSmiles(smi)
    opts = StereoEnumerationOptions(onlyUnassigned=True, unique=True, maxIsomers=limit)
    seen, out = set(), []
    for iso in EnumerateStereoisomers(base, options=opts):
        s = Chem.MolToSmiles(iso)
        mirror = Chem.Mol(iso)
        for a in mirror.GetAtoms():
            if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                a.InvertChirality()
        if Chem.MolToSmiles(mirror) in seen:
            continue
        seen.add(s); out.append(s)
    return out


def model(smi, nconf):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    ps = AllChem.ETKDGv3(); ps.randomSeed = 42; ps.pruneRmsThresh = 0.25; ps.maxIterations = 300
    cids = list(AllChem.EmbedMultipleConfs(m, nconf, ps))
    if not cids:
        return None
    res = AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=20000, mmffVariant='MMFF94s')
    E = np.array([e for _, e in res])
    k = int(np.argmin(E))
    chk = Chem.Mol(m); Chem.AssignStereochemistryFrom3D(chk, confId=cids[k])
    if Chem.MolToSmiles(Chem.RemoveHs(chk)) != Chem.MolToSmiles(Chem.RemoveHs(m)):
        return None                              # optimisation inverted a centre: reject
    w = np.exp(-(E - E.min()) / 0.593); w /= w.sum()
    return m, cids, w, float(E.min())


def score(smi):
    cfg = CFG
    r = model(smi, cfg['nconf'])
    if r is None:
        return None
    m, cids, w, E = r
    idx = {a.GetAtomMapNum(): a.GetIdx() for a in m.GetAtoms() if a.GetAtomMapNum()}
    H = lambda c: [n.GetIdx() for n in m.GetAtomWithIdx(idx[c]).GetNeighbors() if n.GetAtomicNum() == 1]
    confs = [m.GetConformer(c) for c in cids]

    def J(a, b, ha, hb):
        return float(np.dot(w, [karplus.j_hla(m, cf, ha, idx[a], idx[b], hb) for cf in confs]))

    def r6(ha, hb):
        v = np.array([T.GetBondLength(cf, ha, hb) for cf in confs])
        return float(np.dot(w, v ** -6.0) ** (-1 / 6))
    s, detail = 0.0, {}
    for a, b, jexp in cfg.get('j', []):
        vals = [J(a, b, ha, hb) for ha in H(a) for hb in H(b)]
        jc = min(vals, key=lambda v: abs(v - jexp))
        s += (jc - jexp) ** 2; detail['J%d,%d' % (a, b)] = round(jc, 1)
    for a, bs, jexp in cfg.get('j_sets', []):
        calc = sorted([J(a, b, H(a)[0], hb) for b in bs for hb in H(b)], reverse=True)
        expv = sorted(jexp, reverse=True) + [0.0] * (len(calc) - len(jexp))
        s += sum((c - e) ** 2 for c, e in zip(calc, expv))
        detail['J%d,{%s}' % (a, ','.join(map(str, bs)))] = '/'.join('%.1f' % c for c in calc)
    wn, rmax, rmin = cfg.get('w_noe', 4), cfg.get('r_max', 3.2), cfg.get('r_min', 3.0)
    for a, b in cfg.get('noe', []):
        d = min(r6(ha, hb) for ha in H(a) for hb in H(b))
        s += wn * max(0, d - rmax) ** 2; detail['r%d-%d' % (a, b)] = round(d, 1)
    for a, b in cfg.get('no_noe', []):
        d = min(r6(ha, hb) for ha in H(a) for hb in H(b))
        s += wn * max(0, rmin - d) ** 2; detail['r%d-%d(absent)' % (a, b)] = round(d, 1)
    return s, E, smi, detail


def init(cfg):
    global CFG
    CFG = cfg


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('config'); ap.add_argument('--nconf', type=int, default=40)
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 4); ap.add_argument('--top', type=int, default=10)
    a = ap.parse_args()
    cfg = json.load(open(a.config, encoding='utf-8')); cfg['nconf'] = a.nconf
    iso = isomers(cfg['smiles'], cfg.get('max_isomers', 512))
    print('# %d stereoisomers (enantiomers collapsed)' % len(iso), flush=True)
    with Pool(a.workers, initializer=init, initargs=(cfg,)) as p:
        res = [r for r in p.map(score, iso, chunksize=1) if r]
    print('# %d embeddable; %d geometrically impossible' % (len(res), len(iso) - len(res)))
    res.sort(key=lambda r: r[0])
    for s, E, smi, det in res[:a.top]:
        print('score %7.2f  E %7.1f  %s' % (s, E, smi))
        print('              ' + '  '.join('%s=%s' % kv for kv in det.items()))
    if len(res) > 1:
        print('# gap between best and second: %.2f (judge robustness by re-weighting w_noe and '
              'by dropping each observable once)' % (res[1][0] - res[0][0]))


if __name__ == '__main__':
    main()
