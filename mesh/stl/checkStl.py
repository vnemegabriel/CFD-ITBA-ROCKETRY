#!/usr/bin/env python3
"""Checks a rocket STL before meshing and derives the meshInfo reference values.

    ./checkStl.py rocket.stl                       report
    ./checkStl.py rocket.stl --meshInfo meshInfo   report and write the values

Required: solids nosecone, body, boattail and fins; metres; axis along x with
the nose tip at x = 0; four fins in the planes y = 0 and z = 0; a closed
surface.  Exits non-zero, writing nothing, if any of that fails.
"""
import argparse, math, re, sys
from collections import Counter, defaultdict

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('stl')
p.add_argument('--meshInfo')
a = p.parse_args()

solids, cur, tri = defaultdict(list), None, []
with open(a.stl, errors='replace') as f:
    if not f.readline().lstrip().startswith('solid'):
        sys.exit(f'{a.stl}: not an ASCII STL')
    f.seek(0)
    for line in f:
        s = line.split()
        if not s: continue
        if s[0] == 'solid': cur = s[1] if len(s) > 1 else ''
        elif s[0] == 'vertex': tri.append(tuple(map(float, s[1:4])))
        elif s[0] == 'endfacet':
            solids[cur].append(tuple(tri)); tri = []

errors, warnings = [], []
need = ('nosecone', 'body', 'boattail', 'fins')
missing = [n for n in need if n not in solids]
extra = [n for n in solids if n not in need]
if missing or extra:
    sys.exit(f'refused: solids are {list(solids)}; need exactly {list(need)}')

pts = lambda n: [v for t in solids[n] for v in t]
r = lambda v: math.hypot(v[1], v[2])
allp = [v for n in need for v in pts(n)]
xmin, xmax = min(v[0] for v in allp), max(v[0] for v in allp)
rBody = max(r(v) for v in pts('body'))
L = xmax - xmin

if not 0.3 < L < 20:
    errors.append(f'length {L:g}: not metres? (mm would read ~1000x)')
if abs(xmin) > 1e-3 * L:
    errors.append(f'nose tip at x = {xmin:g}, expected 0')
if min(r(v) for v in pts('nosecone')) > 1e-3 * rBody:
    errors.append('the nose has no point on the x axis: the axis is not x')
for k, name in ((1, 'y'), (2, 'z')):
    c = 0.5 * (max(v[k] for v in pts('body')) + min(v[k] for v in pts('body')))
    if abs(c) > 1e-3 * rBody:
        errors.append(f'body off axis in {name} by {c:g}')
off = max(min(abs(v[1]), abs(v[2])) for v in pts('fins'))
if off > 0.25 * rBody:
    errors.append(f'fins not in the planes y = 0 and z = 0 (a fin point sits {off:g} off both)')

key = lambda v: tuple(round(c, 9) for c in v)
edges = Counter()
degenerate = 0
for n in need:
    for t in solids[n]:
        u = [t[1][i] - t[0][i] for i in range(3)]
        w = [t[2][i] - t[0][i] for i in range(3)]
        c = (u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0])
        if math.sqrt(sum(x * x for x in c)) < 1e-14:
            degenerate += 1
        for i in range(3):
            edges[frozenset((key(t[i]), key(t[(i + 1) % 3])))] += 1
bad = sum(1 for n in edges.values() if n != 2)
if bad:
    errors.append(f'surface not closed: {bad} edges not shared by exactly two triangles')
if degenerate:
    errors.append(f'{degenerate} degenerate triangles')

vol = sum(t[0][0] * (t[1][1] * t[2][2] - t[1][2] * t[2][1])
          + t[0][1] * (t[1][2] * t[2][0] - t[1][0] * t[2][2])
          + t[0][2] * (t[1][0] * t[2][1] - t[1][1] * t[2][0])
          for n in need for t in solids[n]) / 6
if vol <= 0:
    warnings.append('normals mostly point inward (cfMesh copes, but check the export)')

fins = pts('fins')
info = {'Aref': math.pi * rBody ** 2, 'lRef': 2 * rBody, 'rBody': rBody, 'lBody': max(v[0] for v in pts('boattail')),
        'finRootR': min(r(v) for v in fins), 'finTipR': max(max(abs(v[1]), abs(v[2])) for v in fins),
        'finX0': min(v[0] for v in fins), 'finX1': max(v[0] for v in fins)}

print(f'{a.stl}: {sum(map(len, solids.values()))} triangles, L = {L:.4f} m, D = {2 * rBody:.4f} m')
for k, v in info.items():
    print(f'  {k:10s} {v:.6g}')
for w in warnings:
    print('  warning:', w)
if errors:
    sys.exit('refused:\n  ' + '\n  '.join(errors))

if a.meshInfo:
    txt = open(a.meshInfo).read()
    for k, v in info.items():
        txt, n = re.subn(rf'^({k}\s+)[^;]*;', rf'\g<1>{v:.6g};', txt, flags=re.M)
        if n != 1:
            sys.exit(f'{a.meshInfo}: no single entry {k}')
    open(a.meshInfo, 'w').write(txt)
    print(f'  -> {a.meshInfo}')
