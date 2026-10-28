"""Haasnoot-de Leeuw-Altona (HLA) generalised Karplus equation for 3J(H,H).

J = P1 cos^2(phi) + P2 cos(phi) + P3
    + sum_i dchi_i * (P4 + P5 * cos^2(xi_i * phi + P6 * |dchi_i|))

dchi_i is the Huggins electronegativity difference (substituent - H), corrected for
beta substituents: dchi_i = dchi_alpha - P7 * sum(dchi_beta).
xi_i = +1/-1 encodes the orientation of substituent i relative to the coupled proton
on the same carbon (see `xi_sign`). The sign convention used here was selected by the
validation in tests/test_karplus.py (pyranoside couplings), not assumed.
Parameter set: 7-parameter general set (Haasnoot et al., Tetrahedron 1980, 36, 2783).
"""
import math
from rdkit.Chem import rdMolTransforms as T

P = dict(P1=13.86, P2=-0.81, P3=0.0, P4=0.56, P5=-2.32, P6=17.9, P7=0.14)
CHI = {'H': 2.20, 'C': 2.60, 'N': 3.05, 'O': 3.50, 'S': 2.58, 'F': 3.90, 'Cl': 3.15, 'Br': 2.96}
XI_CONVENTION = +1   # set by tests/test_karplus.py


def _dchi(mol, sub_idx, carbon_idx):
    """Electronegativity of a substituent, corrected for its own beta substituents."""
    atom = mol.GetAtomWithIdx(sub_idx)
    d = CHI.get(atom.GetSymbol(), 2.6) - CHI['H']
    beta = 0.0
    for nb in atom.GetNeighbors():
        if nb.GetIdx() == carbon_idx or nb.GetAtomicNum() == 1:
            continue
        beta += CHI.get(nb.GetSymbol(), 2.6) - CHI['H']
    return d - P['P7'] * beta


def _proj_angle(conf, x, c1, c2, ref):
    """Signed angle (deg) of atom x around the c1->c2 axis, measured from atom `ref`."""
    import numpy as np
    p = lambda i: np.array(conf.GetAtomPosition(i))
    axis = p(c2) - p(c1); axis /= np.linalg.norm(axis)
    def perp(i):
        v = p(i) - p(c1); return v - axis * np.dot(v, axis)
    u, w = perp(ref), perp(x)
    ang = math.degrees(math.atan2(np.dot(np.cross(u, w), axis), np.dot(u, w)))
    return ang


def j_hla(mol, conf, h1, c1, c2, h2, convention=None):
    """3J(h1,h2) for fragment h1-c1-c2-h2 in conformer `conf` (RDKit, explicit H)."""
    conv = XI_CONVENTION if convention is None else convention
    phi = T.GetDihedralDeg(conf, h1, c1, c2, h2)
    rad = math.radians(phi)
    J = P['P1'] * math.cos(rad) ** 2 + P['P2'] * math.cos(rad) + P['P3']
    for c_same, c_other, h_same in ((c1, c2, h1), (c2, c1, h2)):
        for nb in mol.GetAtomWithIdx(c_same).GetNeighbors():
            s = nb.GetIdx()
            if s in (c_other, h_same) or nb.GetAtomicNum() == 1:
                continue
            dchi = _dchi(mol, s, c_same)
            # orientation of substituent relative to the proton on the same carbon,
            # viewed along c_same -> c_other: positive angle = counter-clockwise (right-hand rule)
            xi = conv * (1 if _proj_angle(conf, s, c_same, c_other, h_same) > 0 else -1)
            J += dchi * (P['P4'] + P['P5'] * math.cos(math.radians(xi * phi + P['P6'] * abs(dchi))) ** 2)
    return J
