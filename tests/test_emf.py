"""EMF parser and axis calibration on synthetic data (no external files).

Run: python -I tests/test_emf.py
"""
import os, struct, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
import emf, spectra


def rec(t, body):
    data = struct.pack('<II', t, 8 + len(body)) + body
    return data + b'\0' * (-len(data) % 4)


def synthetic_emf():
    pen = struct.pack('<5I', 1, 0, 0, 0, 0) + struct.pack('<6I', 0, 1, 0, 0x0000FF, 0, 0)
    poly = struct.pack('<4i', 0, 0, 10, 10) + struct.pack('<I', 3) + struct.pack('<6h', 0, 0, 5, 10, 10, 0)
    s = 'abc'.encode('utf-16-le')
    text = (struct.pack('<4i', 0, 0, 0, 0) + struct.pack('<I2f', 1, 1.0, 1.0)
            + struct.pack('<2i', 7, 9) + struct.pack('<3I', 3, 76, 0) + struct.pack('<4i', 0, 0, 0, 0)
            + struct.pack('<I', 0) + s)
    return rec(1, b'\0' * 80) + rec(95, pen) + rec(37, struct.pack('<I', 1)) + rec(87, poly) + rec(84, text) + rec(14, b'\0' * 12)


def synthetic_spectrum(peak_ppm=3.0):
    """x axis 10..0 ppm over X 0..1000, y axis 0..100 over Y 490..10, one Lorentzian."""
    Y0, X0 = 500.0, 1010.0
    lines = [{'t': 87, 'c': 0, 'p': [[0, Y0], [1000, Y0]]}, {'t': 87, 'c': 0, 'p': [[X0, 0], [X0, 490]]}]
    texts = []
    for k in range(11):
        x = 100.0 * k
        lines.append({'t': 87, 'c': 0, 'p': [[x, Y0], [x, Y0 + 5]]})
        texts.append({'s': '%d' % (10 - k), 'p': [x - 10, Y0 + 20]})
    for k in range(6):
        y = 490 - 96.0 * k
        lines.append({'t': 87, 'c': 0, 'p': [[X0, y], [X0 + 5, y]]})
        texts.append({'s': '%d' % (20 * k), 'p': [X0 + 13, y + 2]})
    X = np.linspace(0, 1000, 2001); ppm = 10 - X / 100
    I = 90 / (1 + ((ppm - peak_ppm) / 0.02) ** 2)
    lines.append({'t': 85, 'c': 0x000080, 'p': [[float(x), float(490 - 4.8 * i)] for x, i in zip(X, I)]})
    return {'lines': lines, 'texts': texts}


def main():
    d = emf.parse_emf(synthetic_emf())
    assert d['lines'][0]['c'] == 0x0000FF, d['lines'][0]
    assert d['lines'][0]['p'] == [(0, 0), (5, 10), (10, 0)], d['lines'][0]['p']
    assert d['texts'][0]['s'] == 'abc' and d['texts'][0]['p'] == (7, 9), d['texts'][0]
    print('emf parse: ok')

    s = synthetic_spectrum(3.0)
    (a, rx), (b, ry), box = spectra.axes(s)
    assert rx < 1e-6 and ry < 1e-6, (rx, ry)
    x, hi, _ = spectra.trace1d(s)
    top = x[int(np.argmax(hi))]
    assert abs(top - 3.0) < 0.01, top
    assert abs(hi.max() - 90) < 0.5, hi.max()
    print('axis calibration and 1D peak: ok (peak %.3f ppm, height %.1f)' % (top, hi.max()))
    print('PASS')


if __name__ == '__main__':
    main()
