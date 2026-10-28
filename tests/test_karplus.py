"""Validate the HLA Karplus implementation and fix its xi sign convention.

Reference values: well-established 3J(H,H) of D-glucopyranose and D-mannopyranose
anomers in D2O (textbook values, +/-0.3 Hz). Run: python -I tests/test_karplus.py
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
from rdkit import Chem
from rdkit.Chem import AllChem
import karplus

# atom maps: ring carbons C1..C5 = maps 1..5
SUGARS = {
    'alpha-glucose':  ('O[C@H:1]1O[C@H:5](CO)[C@@H:4](O)[C@H:3](O)[C@H:2]1O',  {(1, 2): 3.8, (2, 3): 9.6, (3, 4): 9.3, (4, 5): 9.6}),
    'beta-glucose':   ('O[C@@H:1]1O[C@H:5](CO)[C@@H:4](O)[C@H:3](O)[C@H:2]1O', {(1, 2): 8.0, (2, 3): 9.4, (3, 4): 9.2, (4, 5): 9.6}),
    'alpha-mannose':  ('O[C@H:1]1O[C@H:5](CO)[C@@H:4](O)[C@H:3](O)[C@@H:2]1O', {(1, 2): 1.8, (2, 3): 3.4, (3, 4): 9.6}),
    'beta-mannose':   ('O[C@@H:1]1O[C@H:5](CO)[C@@H:4](O)[C@H:3](O)[C@@H:2]1O', {(1, 2): 0.9, (2, 3): 3.4, (3, 4): 9.5}),
}


def predict(smi, pairs, convention, nconf=30):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    ps = AllChem.ETKDGv3(); ps.randomSeed = 7
    cids = list(AllChem.EmbedMultipleConfs(m, nconf, ps))
    res = AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=5000)
    E = np.array([e for _, e in res]); w = np.exp(-(E - E.min()) / 0.593); w /= w.sum()
    idx = {a.GetAtomMapNum(): a.GetIdx() for a in m.GetAtoms() if a.GetAtomMapNum()}
    H = lambda c: [n.GetIdx() for n in m.GetAtomWithIdx(idx[c]).GetNeighbors() if n.GetAtomicNum() == 1][0]
    out = {}
    for a, b in pairs:
        out[(a, b)] = float(np.dot(w, [karplus.j_hla(m, m.GetConformer(c), H(a), idx[a], idx[b], H(b), convention) for c in cids]))
    return out


def rmsd(convention, verbose=False):
    err = []
    for name, (smi, ref) in SUGARS.items():
        pr = predict(smi, ref, convention)
        for k, v in ref.items():
            err.append(pr[k] - v)
            if verbose:
                print('  %-14s J%d,%d  calc %5.1f  exp %4.1f' % (name, *k, pr[k], v))
    return float(np.sqrt(np.mean(np.square(err))))


if __name__ == '__main__':
    scores = {c: rmsd(c) for c in (+1, -1)}
    for c, s in scores.items():
        print('xi convention %+d: RMSD %.2f Hz' % (c, s))
    best = min(scores, key=scores.get)
    print('module XI_CONVENTION = %+d ; best = %+d' % (karplus.XI_CONVENTION, best))
    rmsd(best, verbose=True)
    assert best == karplus.XI_CONVENTION, 'XI_CONVENTION in karplus.py does not match validation'
    assert scores[best] < 1.0, 'RMSD above 1 Hz: implementation problem'
    print('PASS')
