"""Test candidate constitutions against observed HMBC / COSY correlations.

    python topology.py CONFIG.json

For every candidate (atom-mapped SMILES) all assignments of observed signals to atoms of
the same class are tried (permutations inside each group). Each assignment is scored by the
number of violated correlations:
  HMBC  violated if the shortest path H-bearing carbon -> target carbon is > max_hmbc_bonds - 1
        bonds (default max 3J, i.e. carbons at most 2 bonds apart)
  COSY  violated if the two H-bearing carbons are not directly bonded (3J H,H)
Report: best violation count, how many assignments reach it, the violations and the mapping.
A candidate that cannot reach zero needs every remaining violation explained explicitly.

CONFIG keys
  candidates : {name: mapped SMILES}
  groups     : {group: {"signals": [ids], "atoms": [maps] or {candidate: [maps]}}}
  fixed      : {signal: map or {candidate: map}}
  aliases    : {id: [signal, ...]}   overlapped 1H signal that may sit on any listed carbon
  hmbc       : [[H_signal, [C_signal, ...alternatives]], ...]
  cosy       : [[H_signal, H_signal], ...]
  max_hmbc_bonds (3), max_assignments (2e6)
Signal ids are free strings (e.g. "c1", or the carbon shift as text); they only have to be consistent.
"""
import itertools, json, math, sys
import numpy as np
from rdkit import Chem


def _per_cand(v, name):
    return v[name] if isinstance(v, dict) else v


def evaluate(cfg, name, smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError('bad SMILES for %s' % name)
    idx = {a.GetAtomMapNum(): a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum()}
    D = Chem.GetDistanceMatrix(mol)
    sig, perms = [], []
    for g, spec in cfg['groups'].items():
        atoms = _per_cand(spec['atoms'], name)
        if len(atoms) != len(spec['signals']):
            raise ValueError('group %s: %d signals vs %d atoms for %s' % (g, len(spec['signals']), len(atoms), name))
        sig += spec['signals']
        perms.append(list(itertools.permutations(atoms)))
    fixed = {s: _per_cand(m, name) for s, m in cfg.get('fixed', {}).items()}
    sig += list(fixed)
    n = math.prod(len(p) for p in perms)
    if n > cfg.get('max_assignments', 2e6):
        raise ValueError('%s: %d assignments; split groups or fix more signals' % (name, n))
    A = np.array([[idx[m] for p in combo for m in p] + [idx[m] for m in fixed.values()]
                  for combo in itertools.product(*perms)])
    col = {s: i for i, s in enumerate(sig)}
    alias = cfg.get('aliases', {})

    def atoms_of(s):
        return [A[:, col[x]] for x in alias.get(s, [s])]
    viol = np.zeros(len(A), int)
    checks = []
    maxd = cfg.get('max_hmbc_bonds', 3) - 1
    for h, targets in cfg.get('hmbc', []):
        targets = targets if isinstance(targets, list) else [targets]
        d = np.min([D[a, A[:, col[t]]] for a in atoms_of(h) for t in targets], axis=0)
        viol += d > maxd; checks.append(('HMBC', h, targets, d))
    for x, y in cfg.get('cosy', []):
        d = np.min([D[a, b] for a in atoms_of(x) for b in atoms_of(y)], axis=0)
        viol += d != 1; checks.append(('COSY', x, y, d))
    k = int(np.argmin(viol)); best = int(viol[k])
    print('== %s: %d assignments, fewest violations %d (reached by %d)'
          % (name, len(A), best, int((viol == best).sum())))
    inv = {v: k2 for k2, v in idx.items()}
    for kind, a, b, d in checks:
        if (kind == 'HMBC' and d[k] > maxd) or (kind == 'COSY' and d[k] != 1):
            nJ = int(d[k]) + (1 if kind == 'HMBC' else 2)
            print('   violated %s %s -> %s : %dJ' % (kind, a, b, nJ))
    print('   mapping:', ', '.join('%s=C%d' % (s, inv[A[k, col[s]]]) for s in sig))
    return best


if __name__ == '__main__':
    cfg = json.load(open(sys.argv[1], encoding='utf-8'))
    res = {n: evaluate(cfg, n, s) for n, s in cfg['candidates'].items()}
    print('\nsummary:', ', '.join('%s=%d' % kv for kv in sorted(res.items(), key=lambda kv: kv[1])))
