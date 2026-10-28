"""Calibrate axes and extract peaks from EMF-derived spectrum JSON (output of emf.py).

Axis detection is automatic: the longest horizontal and vertical black (pen colour 0)
segments are taken as the x and y axes; tick marks are short black segments that start on
an axis; numeric text labels next to the axis are matched to ticks with a fitted offset.
The fit residual is printed so you can judge the calibration. Nothing about one specific
export layout is hard-coded.

    python spectra.py peaks1d  FILE.json [--thr 0.02] [--both]     1D peak list (ppm, height)
    python spectra.py integ    FILE.json --ref LO HI NH --reg LO HI [--reg LO HI ...]
    python spectra.py mult     FILE.json LO HI                     multiplet lines + gaps (Hz)
    python spectra.py peaks2d  FILE.json [--gap 0.5] [--out F.tsv] 2D cross-peaks (F2, F1)
    python spectra.py proj     FILE.json [--thr 0.05]              1D projections of a 2D plot
    python spectra.py zoom     FILE.json F2a F2b F1a F1b OUT.png   re-plot a 2D region

Conventions: ppm axes decrease to the right / downward as drawn. Spectrometer frequency for
Hz conversion is --mhz (default 400.13). Colours are EMF COLORREF ints (0x00BBGGRR).
"""
import argparse, json, sys
import numpy as np

AXIS = 0                 # black
GRID = 0xC0C0C0          # light grey grid
TRACE = 0x000080         # dark red 1D traces (MestReNova default); auto-detected if absent
MARKER = 0xFFFF00        # cyan cursor/marker lines drawn by MestReNova; never a contour


def load(fn):
    with open(fn) as f:
        return json.load(f)


def _segs(d, colour=AXIS):
    return [l['p'] for l in d['lines'] if l['c'] == colour and len(l['p']) == 2]


def axes(d):
    """Return ((a, res_x), (b, res_y), box). ppm_x = polyval(a, X); value_y = polyval(b, Y)."""
    L = _segs(d)
    hz = [p for p in L if abs(p[0][1] - p[1][1]) < .01]
    vt = [p for p in L if abs(p[0][0] - p[1][0]) < .01]
    hx = max(hz, key=lambda p: abs(p[1][0] - p[0][0]))
    vy = max(vt, key=lambda p: abs(p[1][1] - p[0][1]))
    Y0, X0 = hx[0][1], vy[0][0]
    span_x = abs(hx[1][0] - hx[0][0]); span_y = abs(vy[1][1] - vy[0][1])
    def major(ticks):
        """Keep only major ticks (the longest tick class); minor ticks would let labels
        fit equally well at half-interval offsets."""
        if not ticks:
            return []
        Lmax = max(L for _, L in ticks)
        return [pos for pos, L in ticks if L > .75 * Lmax]
    xt = major([(p[0][0], abs(p[1][1] - p[0][1])) for p in vt
                if abs(min(p[0][1], p[1][1]) - Y0) < .6 and abs(p[1][1] - p[0][1]) < .05 * span_y])
    yt = major([(p[0][1], abs(p[1][0] - p[0][0])) for p in hz
                if abs(min(p[0][0], p[1][0]) - X0) < .6 and abs(p[1][0] - p[0][0]) < .05 * span_x])

    def labels(axis):
        """Numeric labels in the single row (x) / column (y) closest to the axis. Rotated
        axis titles such as 'f1 (ppm)' are drawn char by char; '1' would parse as a number,
        so only the most populated row/column is used."""
        from collections import Counter
        cand = []
        for t in d['texts']:
            try:
                v = float(t['s'])
            except ValueError:
                continue
            if axis == 'x' and Y0 < t['p'][1] < Y0 + .08 * span_y:
                cand.append((round(t['p'][1]), t['p'][0], v))
            if axis == 'y' and X0 < t['p'][0] < X0 + .08 * span_x:
                cand.append((round(t['p'][0]), t['p'][1], v))
        if not cand:
            return []
        row = Counter(c[0] for c in cand).most_common(1)[0][0]
        return [(pos, v) for r, pos, v in cand if abs(r - row) <= 1]

    def fit(ticks, labs, span):
        best = None
        for off in np.linspace(-.03 * span, .03 * span, 121):
            pairs = [(min(ticks, key=lambda q: abs(q - (p + off))), v) for p, v in labs]
            if len({p for p, _ in pairs}) < max(3, len(pairs) - 1):
                continue
            a = np.polyfit(*zip(*pairs), 1)
            r = max(abs(np.polyval(a, p) - v) for p, v in pairs)
            if best is None or r < best[1]:
                best = (a, r)
        return best

    fx = fit(xt, labels('x'), span_x)
    fy = fit(yt, labels('y'), span_y) if yt and labels('y') else (None, None)
    box = (min(hx[0][0], hx[1][0]), max(hx[0][0], hx[1][0]), min(vy[0][1], vy[1][1]), max(vy[0][1], vy[1][1]))
    return fx, fy, box


def trace_colour(d):
    """Colour of the long 1D trace polylines (most points among non-axis colours)."""
    from collections import Counter
    c = Counter()
    for l in d['lines']:
        if l['c'] not in (AXIS, GRID, None) and len(l['p']) > 50:
            c[l['c']] += len(l['p'])
    return c.most_common(1)[0][0] if c else TRACE


def trace1d(d):
    """(ppm, intensity) of the main 1D trace, one max/min envelope point per unique x."""
    (a, _), (b, _), box = axes(d)
    col = trace_colour(d)
    P = np.vstack([np.array(l['p']) for l in d['lines'] if l['c'] == col and len(l['p']) > 50])
    x = np.polyval(a, P[:, 0]); y = np.polyval(b, P[:, 1]) if b is not None else -P[:, 1]
    xs, inv = np.unique(np.round(x, 5), return_inverse=True)
    hi = np.full(len(xs), -np.inf); lo = np.full(len(xs), np.inf)
    np.maximum.at(hi, inv, y); np.minimum.at(lo, inv, y)
    return xs, hi, lo


def peaks1d(fn, thr=0.02, both=False):
    d = load(fn); (a, rx), (b, ry), _ = axes(d)
    x, hi, lo = trace1d(d)
    print('# calibration residual: x %.4f ppm, y %s ; points %d ; median step %.4f ppm'
          % (rx, 'n/a' if ry is None else '%.3g' % ry, len(x), np.median(np.diff(x))))
    top = max(abs(hi).max(), abs(lo).max())
    out = []
    for sgn, arr in ((1, hi), (-1, lo)) if both else ((1, hi),):
        yy = sgn * arr
        for i in range(2, len(x) - 2):
            if yy[i] > thr * top and yy[i] == yy[max(0, i - 3):i + 4].max():
                out.append((x[i], sgn * yy[i]))
    for p, h in sorted(out, reverse=True):
        print('%9.3f  %+12.4g' % (p, h))


def integ(fn, ref, regs):
    d = load(fn); x, hi, _ = trace1d(d)
    I = lambda h, l: np.trapezoid(hi[(x <= max(h, l)) & (x >= min(h, l))], x[(x <= max(h, l)) & (x >= min(h, l))])
    unit = I(ref[0], ref[1]) / ref[2]
    print('# reference %.2f-%.2f = %g H. Note: decimated traces under-integrate sharp singlets.' % ref)
    for h, l in regs:
        print('%6.2f-%6.2f  %6.2f H' % (h, l, I(h, l) / unit))


def mult(fn, lo, hi, mhz=400.13, rel=0.15):
    d = load(fn); x, y, _ = trace1d(d)
    m = (x >= min(lo, hi)) & (x <= max(lo, hi)); X, Y = x[m], y[m]
    pk = [X[i] for i in range(1, len(X) - 1) if Y[i] > Y[i - 1] and Y[i] >= Y[i + 1] and Y[i] > rel * Y.max()]
    pk = sorted(pk, reverse=True)
    print('lines (ppm):', ' '.join('%.4f' % p for p in pk))
    print('gaps (Hz):  ', ' '.join('%.1f' % ((pk[i] - pk[i + 1]) * mhz) for i in range(len(pk) - 1)))
    print('# dd: gaps = J2, J1-J2, J2. Digital resolution limits J to about +/-0.5 Hz.')


def contours(d, box):
    skip = {AXIS, GRID, MARKER, None, trace_colour(d)}
    out = []
    for l in d['lines']:
        if l['c'] in skip or len(l['p']) < 2:
            continue
        p = np.array(l['p'])
        if p[:, 0].min() < box[0] - 1 or p[:, 0].max() > box[1] + 1 or p[:, 1].min() < box[2] - 1 or p[:, 1].max() > box[3] + 1:
            continue
        out.append((l['c'], p))
    return out


def cluster(cs, gap):
    """Union-find on overlapping bounding boxes (dilated by `gap` device units)."""
    bb = np.array([[p[:, 0].min(), p[:, 0].max(), p[:, 1].min(), p[:, 1].max()] for _, p in cs])
    par = list(range(len(cs)))

    def f(i):
        while par[i] != i:
            par[i] = par[par[i]]; i = par[i]
        return i
    o = np.argsort(bb[:, 0])
    for k, i in enumerate(o):
        for j in o[k + 1:]:
            if bb[j, 0] > bb[i, 1] + gap:
                break
            if bb[j, 2] <= bb[i, 3] + gap and bb[i, 2] <= bb[j, 3] + gap:
                par[f(i)] = f(j)
    g = {}
    for i in range(len(cs)):
        g.setdefault(f(i), []).append(i)
    return list(g.values())


def peaks2d(fn, gap=0.5, out=None):
    """Cross-peak list. Centre = centre of the innermost (smallest) contour colour level.
    Colours tell the sign: in MestReNova defaults positive levels are red/orange/yellow and
    negative levels blue. Report them; interpret per experiment (e.g. NOESY vs diagonal)."""
    d = load(fn); (a, rx), (b, ry), box = axes(d)
    cs = contours(d, box)
    rows = []
    for g in cluster(cs, gap):
        bycol = {}
        for i in g:
            bycol.setdefault(cs[i][0], []).append(i)
        P = min((np.vstack([cs[i][1] for i in ii]) for ii in bycol.values()),
                key=lambda P: (np.ptp(P[:, 0]) + .5) * (np.ptp(P[:, 1]) + .5))
        cx, cy = (P[:, 0].min() + P[:, 0].max()) / 2, (P[:, 1].min() + P[:, 1].max()) / 2
        Q = np.vstack([cs[i][1] for i in g])
        rows.append((np.polyval(a, cx), np.polyval(b, cy), len(g),
                     np.polyval(a, Q[:, 0].max()), np.polyval(a, Q[:, 0].min()),
                     np.polyval(b, Q[:, 1].min()), np.polyval(b, Q[:, 1].max()),
                     ','.join('%06X' % c for c in sorted(bycol))))
    rows.sort(key=lambda r: (-r[0], r[1]))
    lines = ['# calibration residual F2 %.4f ppm, F1 %.4f ppm; n_contours = strength proxy' % (rx, ry),
             '# F2\tF1\tn\tF2_lo\tF2_hi\tF1_lo\tF1_hi\tcolours']
    lines += ['%.3f\t%.2f\t%d\t%.3f\t%.3f\t%.2f\t%.2f\t%s' % r for r in rows]
    txt = '\n'.join(lines)
    if out:
        with open(out, 'w') as f:
            f.write(txt + '\n')
    print(txt)


def proj(fn, thr=0.05):
    """Peak positions of the 1D traces drawn along the top and left of a 2D plot. Use them to
    check that 1D and 2D data come from the same sample, and as an extra 13C peak list."""
    d = load(fn); (a, _), (b, _), box = axes(d)
    col = trace_colour(d)
    for side in ('top', 'left'):
        P = [np.array(l['p']) for l in d['lines'] if l['c'] == col and len(l['p']) > 50]
        P = [p for p in P if (p[:, 1].max() < box[2] if side == 'top' else p[:, 0].max() < box[0])]
        if not P:
            print(side, ': none'); continue
        P = np.vstack(P)
        if side == 'top':
            ppm, I = np.polyval(a, P[:, 0]), np.median(P[:, 1]) - P[:, 1]
        else:
            ppm, I = np.polyval(b, P[:, 1]), np.median(P[:, 0]) - P[:, 0]
        u, inv = np.unique(np.round(ppm, 4), return_inverse=True)
        m = np.full(len(u), -np.inf); np.maximum.at(m, inv, I)
        pk = [u[i] for i in range(2, len(u) - 2) if m[i] > thr * m.max() and m[i] == m[max(0, i - 3):i + 4].max()]
        print('%s (step %.3f ppm): %s' % (side, np.median(np.diff(u)), ' '.join('%.2f' % v for v in pk)))


def zoom(fn, f2a, f2b, f1a, f1b, png):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    d = load(fn); (a, _), (b, _), box = axes(d)
    fig, ax = plt.subplots(figsize=(7, 5), dpi=110)
    for c, p in contours(d, box):
        X, Y = np.polyval(a, p[:, 0]), np.polyval(b, p[:, 1])
        if X.max() < min(f2a, f2b) or X.min() > max(f2a, f2b) or Y.max() < min(f1a, f1b) or Y.min() > max(f1a, f1b):
            continue
        ax.plot(X, Y, color='#%02x%02x%02x' % (c & 255, (c >> 8) & 255, (c >> 16) & 255), lw=.6)
    ax.set_xlim(max(f2a, f2b), min(f2a, f2b)); ax.set_ylim(max(f1a, f1b), min(f1a, f1b))
    ax.set_xlabel('F2 (ppm)'); ax.set_ylabel('F1 (ppm)'); ax.grid(alpha=.3); ax.minorticks_on()
    fig.tight_layout(); fig.savefig(png); print('saved', png)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('peaks1d'); p.add_argument('file'); p.add_argument('--thr', type=float, default=.02); p.add_argument('--both', action='store_true')
    p = sp.add_parser('integ'); p.add_argument('file'); p.add_argument('--ref', nargs=3, type=float, required=True)
    p.add_argument('--reg', nargs=2, type=float, action='append', required=True)
    p = sp.add_parser('mult'); p.add_argument('file'); p.add_argument('lo', type=float); p.add_argument('hi', type=float); p.add_argument('--mhz', type=float, default=400.13)
    p = sp.add_parser('peaks2d'); p.add_argument('file'); p.add_argument('--gap', type=float, default=.5); p.add_argument('--out')
    p = sp.add_parser('proj'); p.add_argument('file'); p.add_argument('--thr', type=float, default=.05)
    p = sp.add_parser('zoom'); p.add_argument('file'); [p.add_argument(k, type=float) for k in ('f2a', 'f2b', 'f1a', 'f1b')]; p.add_argument('png')
    a = ap.parse_args()
    if a.cmd == 'peaks1d': peaks1d(a.file, a.thr, a.both)
    elif a.cmd == 'integ': integ(a.file, tuple(a.ref), a.reg)
    elif a.cmd == 'mult': mult(a.file, a.lo, a.hi, a.mhz)
    elif a.cmd == 'peaks2d': peaks2d(a.file, a.gap, a.out)
    elif a.cmd == 'proj': proj(a.file, a.thr)
    elif a.cmd == 'zoom': zoom(a.file, a.f2a, a.f2b, a.f1a, a.f1b, a.png)


if __name__ == '__main__':
    main()
