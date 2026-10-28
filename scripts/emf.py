"""Read vector spectra (EMF) from Office documents or standalone .emf files.

    python emf.py INPUT [--out DIR]

INPUT is a .docx/.pptx file or an .emf file. For every EMF image one JSON file is
written to DIR (default: ./work/emf) with
    lines: [{"t": record type, "c": pen colour (0x00BBGGRR int), "p": [[x, y], ...]}, ...]
    texts: [{"s": string, "p": [x, y]}, ...]
in device coordinates after the EMF world transform. index.json lists the images in
document order with the nearest preceding caption text.

Privacy: embedded OLE objects (e.g. MestReNova documents) are NOT extracted. They often
contain acquisition paths, user names and sample IDs. Never redistribute them or the
source document.
"""
import argparse, json, os, posixpath, re, struct, sys, zipfile


def _records(b):
    off = 0
    while off + 8 <= len(b):
        t, sz = struct.unpack_from('<II', b, off)
        if sz == 0:
            break
        yield t, b[off:off + sz]
        off += sz


def _mul(a, b):
    a11, a12, a21, a22, adx, ady = a
    b11, b12, b21, b22, bdx, bdy = b
    return (a11 * b11 + a12 * b21, a11 * b12 + a12 * b22, a21 * b11 + a22 * b21,
            a21 * b12 + a22 * b22, adx * b11 + ady * b21 + bdx, adx * b12 + ady * b22 + bdy)


def _ap(x, p):
    return (p[0] * x[0] + p[1] * x[2] + x[4], p[0] * x[1] + p[1] * x[3] + x[5])


def parse_emf(data):
    """Return {'lines': [...], 'texts': [...]} from EMF bytes."""
    objs, pen, xf, stack = {}, None, (1, 0, 0, 1, 0, 0), []
    lines, texts = [], []
    for t, r in _records(data):
        if t == 95:                                   # EMR_EXTCREATEPEN
            ih = struct.unpack_from('<I', r, 8)[0]
            color = struct.unpack_from('<I', r, 40)[0]
            objs[ih] = ('pen', color & 0xFFFFFF)
        elif t == 38:                                 # EMR_CREATEPEN
            ih, _s, _wx, _wy, color = struct.unpack_from('<5I', r, 8)
            objs[ih] = ('pen', color & 0xFFFFFF)
        elif t == 37:                                 # EMR_SELECTOBJECT
            ih = struct.unpack_from('<I', r, 8)[0]
            if ih in objs and objs[ih][0] == 'pen':
                pen = objs[ih]
        elif t == 35:                                 # EMR_SETWORLDTRANSFORM
            xf = struct.unpack_from('<6f', r, 8)
        elif t == 36:                                 # EMR_MODIFYWORLDTRANSFORM
            m = struct.unpack_from('<6f', r, 8)
            mode = struct.unpack_from('<I', r, 32)[0]
            xf = {1: (1, 0, 0, 1, 0, 0), 2: _mul(m, xf), 3: _mul(xf, m), 4: m}.get(mode, xf)
        elif t == 33:                                 # EMR_SAVEDC
            stack.append(xf)
        elif t == 34:                                 # EMR_RESTOREDC
            n = struct.unpack_from('<i', r, 8)[0]
            k = -n if n < 0 else len(stack) - n + 1
            for _ in range(max(k, 0)):
                if stack:
                    xf = stack.pop()
        elif t in (85, 86, 87, 88, 89):               # 16-bit poly records
            n = struct.unpack_from('<I', r, 24)[0]
            pts = [_ap(xf, struct.unpack_from('<2h', r, 28 + 4 * i)) for i in range(n)]
            lines.append({'t': t, 'c': pen[1] if pen else None, 'p': pts})
        elif t in (2, 3, 4, 5, 6):                    # 32-bit poly records
            n = struct.unpack_from('<I', r, 24)[0]
            pts = [_ap(xf, struct.unpack_from('<2i', r, 28 + 8 * i)) for i in range(n)]
            lines.append({'t': t, 'c': pen[1] if pen else None, 'p': pts})
        elif t in (90, 7):                            # poly-polyline 16 / 32 bit
            npl = struct.unpack_from('<I', r, 24)[0]
            cnt = struct.unpack_from('<%dI' % npl, r, 32)
            o = 32 + 4 * npl
            w, fmt = (4, '<2h') if t == 90 else (8, '<2i')
            for c in cnt:
                pts = [_ap(xf, struct.unpack_from(fmt, r, o + w * i)) for i in range(c)]
                o += w * c
                lines.append({'t': t, 'c': pen[1] if pen else None, 'p': pts})
        elif t == 84:                                 # EMR_EXTTEXTOUTW
            rx, ry, nch, offs = struct.unpack_from('<2i2I', r, 36)
            s = r[offs:offs + 2 * nch].decode('utf-16-le', 'replace')
            texts.append({'s': s, 'p': _ap(xf, (rx, ry))})
    return {'lines': lines, 'texts': texts}


RE_REL = re.compile(r'<Relationship [^>]*>')
RE_PARA = re.compile(r'<(?:w|a):p[ >].*?</(?:w|a):p>', re.S)
RE_TEXT = re.compile(r'<(?:w|a):t[^>]*>([^<]*)</(?:w|a):t>')


def _office_images(path):
    """Yield (name, bytes, caption) for EMF images in document order."""
    z = zipfile.ZipFile(path)
    names = z.namelist()
    emb = [n for n in names if '/embeddings/' in n]
    if emb:
        print('note: %d embedded OLE object(s) skipped (may contain private paths or user '
              'names; do not redistribute)' % len(emb), file=sys.stderr)
    seen = []
    docs = sorted(n for n in names if re.match(r'(word/document|ppt/slides/slide[0-9]+)[.]xml$', n))
    for dx in docs:
        folder, base = dx.rsplit('/', 1)
        rels = folder + '/_rels/' + base + '.rels'
        if rels not in names:
            continue
        rmap = {}
        for el in RE_REL.findall(z.read(rels).decode('utf-8')):
            i = re.search(r'Id="([^"]+)"', el)
            tg = re.search(r'Target="([^"]+)"', el)
            if i and tg:
                rmap[i.group(1)] = tg.group(1)
        xml = z.read(dx).decode('utf-8')
        caption = ''
        for para in RE_PARA.findall(xml):
            txt = ''.join(RE_TEXT.findall(para)).strip()
            if txt:
                caption = txt
            for rid in re.findall(r'r:(?:embed|id)="([^"]+)"', para):
                tgt = rmap.get(rid, '')
                if tgt.lower().endswith('.emf'):
                    full = posixpath.normpath(posixpath.join(folder, tgt))
                    if full in names and full not in seen:
                        seen.append(full)
                        yield full.rsplit('/', 1)[1], z.read(full), caption
    for n in names:
        if n.lower().endswith('.emf') and n not in seen:
            yield n.rsplit('/', 1)[1], z.read(n), ''


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('input')
    ap.add_argument('--out', default=os.path.join('work', 'emf'))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.input.lower().endswith('.emf'):
        with open(a.input, 'rb') as f:
            items = [(os.path.basename(a.input), f.read(), '')]
    else:
        items = list(_office_images(a.input))
    index = []
    for k, (name, data, cap) in enumerate(items, 1):
        d = parse_emf(data)
        fn = os.path.join(a.out, '%02d_%s.json' % (k, os.path.splitext(name)[0]))
        with open(fn, 'w') as f:
            json.dump(d, f)
        index.append({'order': k, 'file': os.path.basename(fn), 'caption': cap[:120],
                      'n_lines': len(d['lines']), 'n_texts': len(d['texts'])})
        print('%2d  %-22s lines=%-6d texts=%-4d  %s' % (k, os.path.basename(fn), len(d['lines']),
                                                      len(d['texts']), cap[:60]))
    with open(os.path.join(a.out, 'index.json'), 'w', encoding='utf-8') as f:
        json.dump(index, f, indent=1, ensure_ascii=False)


if __name__ == '__main__':
    main()
