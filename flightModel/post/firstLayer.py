#!/usr/bin/env python3
"""firstLayer.py <mesh-dir>

The first cell at each wall patch: y, the distance from the cell centre to
the wall face (what the wall functions and yPlus use), and t1 = 2 y, its
thickness.  Per patch, weighted by face area: the median and the 5 %
quantile, and the minimum.  By area, because a refined band has many small
faces and would pull a plain median towards it.

Runs postProcess -func writeCellCentres on a temporary case linked to the
mesh, so the mesh directory is not touched.  Needs the OpenFOAM environment
and an ASCII mesh.
"""
import math
import os
import re
import subprocess
import sys
import tempfile

if len(sys.argv) != 2:
    sys.exit(__doc__)
mesh = os.path.abspath(sys.argv[1])
pm = os.path.join(mesh, 'constant', 'polyMesh')


def entries(path):
    """The size of the top-level list of an OpenFOAM file, and its lines."""
    f = open(path)
    prev = ''
    for line in f:
        if line.startswith('format') and 'binary' in line:
            sys.exit(f'{path} is binary: the mesh must be written in ascii')
        if line.strip() == '(' and prev.strip().isdigit():
            return int(prev), f
        prev = line
    sys.exit(f'no list in {path}')


def items(path):
    """One entry per line, as many as the list says."""
    n, f = entries(path)
    for _ in range(n):
        yield next(f)


def faces(path):
    """Each face as a list of point labels, in order."""
    n, f = entries(path)
    for _ in range(n):
        s = next(f).strip()
        while not s:
            s = next(f).strip()
        if s.endswith(')'):
            yield [int(v) for v in s[s.index('(') + 1:-1].split()]
        else:
            next(f)
            yield [int(next(f)) for _ in range(int(s))]
            next(f)


boundary = open(os.path.join(pm, 'boundary')).read()
walls = {}
for m in re.finditer(r'\n\s*(\w+)\s*\{([^}]*)\}', boundary):
    body = m.group(2)
    if re.search(r'type\s+wall;', body):
        start = int(re.search(r'startFace\s+(\d+)', body).group(1))
        walls[m.group(1)] = (start, start + int(re.search(r'nFaces\s+(\d+)', body).group(1)))
if not walls:
    sys.exit('no wall patch in ' + pm)
lo, hi = min(a for a, _ in walls.values()), max(b for _, b in walls.values())

wallFaces = {}
for i, f in enumerate(faces(os.path.join(pm, 'faces'))):
    if i >= hi:
        break
    if i >= lo:
        wallFaces[i] = f
owner = {}
for i, line in enumerate(items(os.path.join(pm, 'owner'))):
    if i >= hi:
        break
    if i >= lo:
        owner[i] = int(line)
need = {p for f in wallFaces.values() for p in f}
points = {}
for i, line in enumerate(items(os.path.join(pm, 'points'))):
    if i in need:
        points[i] = tuple(map(float, line.strip().strip('()').split()))

with tempfile.TemporaryDirectory() as case:
    os.makedirs(os.path.join(case, 'constant'))
    os.symlink(pm, os.path.join(case, 'constant', 'polyMesh'))
    os.symlink(os.path.join(mesh, 'system'), os.path.join(case, 'system'))
    subprocess.run(['postProcess', '-case', case, '-func', 'writeCellCentres', '-constant',
                    '-noFunctionObjects'], check=True, stdout=subprocess.DEVNULL)
    cells = set(owner.values())
    centre = {}
    with open(os.path.join(case, 'constant', 'C')) as f:
        for line in f:
            if line.startswith('internalField'):
                break
        n = int(next(f))
        next(f)
        for i in range(n):
            line = next(f)
            if i in cells:
                centre[i] = tuple(map(float, line.strip().strip('()').split()))



def quantile(pairs, q):
    total = sum(w for _, w in pairs)
    acc = 0.0
    for v, w in pairs:
        acc += w
        if acc >= q * total:
            return v
    return pairs[-1][0]


print(f'{mesh}')
print(f'{"patch":8s} {"faces":>8s} {"y median":>10s} {"t1 median":>10s} {"t1 5 %":>9s} {"t1 min":>9s}   [mm, by area]')
for name, (a, b) in walls.items():
    t1 = []
    for i in range(a, b):
        v = [points[p] for p in wallFaces[i]]
        c = [sum(p[d] for p in v) / len(v) for d in range(3)]
        n = [0.0, 0.0, 0.0]
        for p, q in zip(v, v[1:] + v[:1]):
            n[0] += (p[1] - q[1]) * (p[2] + q[2])
            n[1] += (p[2] - q[2]) * (p[0] + q[0])
            n[2] += (p[0] - q[0]) * (p[1] + q[1])
        ln = math.sqrt(sum(x * x for x in n))
        o = centre[owner[i]]
        t1.append((2 * abs(sum((o[d] - c[d]) * n[d] / ln for d in range(3))), ln / 2))
    t1.sort()
    med = quantile(t1, 0.5)
    print(f'{name:8s} {len(t1):8d} {med / 2 * 1e3:10.4f} {med * 1e3:10.4f} '
          f'{quantile(t1, 0.05) * 1e3:9.4f} {t1[0][0] * 1e3:9.4f}')
