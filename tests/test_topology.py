"""topology.py on a toy case: ethyl acetate vs methyl propanoate.

Observed (generic ids): COSY CH2/CH3a; HMBC CH3a -> C=O, CH2 -> C=O, CH3b -> C=O.
CH3a -> C=O is 3J in methyl propanoate but 4J in ethyl acetate.
Run: python -I tests/test_topology.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
import topology

CFG = {
    'candidates': {'ethyl acetate': '[CH3:1][C:2](=O)O[CH2:3][CH3:4]',
                   'methyl propanoate': '[CH3:1][CH2:3][C:2](=O)O[CH3:4]'},
    'groups': {'CH3': {'signals': ['CH3a', 'CH3b'], 'atoms': [1, 4]}},
    'fixed': {'C=O': 2, 'CH2': 3},
    'hmbc': [['CH3a', ['C=O']], ['CH2', ['C=O']], ['CH3b', ['C=O']]],
    'cosy': [['CH2', 'CH3a']],
}

if __name__ == '__main__':
    r = {n: topology.evaluate(CFG, n, s) for n, s in CFG['candidates'].items()}
    assert r == {'ethyl acetate': 1, 'methyl propanoate': 0}, r
    print('PASS')
