#!/usr/bin/env python3
"""snapPlane.py <case> <patch>

Sets z = 0 on every point of <patch> in <case>/constant/polyMesh/points.
The half mesh is the full mesh cut at z = 0 by subsetMesh; the cut faces
lie on the plane within the smoothing of cartesianMesh, and a symmetry patch
must be flat.  Needs writeFormat ascii.
"""
import re
import sys

if len(sys.argv) != 3:
    sys.exit(__doc__)
pm = sys.argv[1] + '/constant/polyMesh/'
patch = sys.argv[2]


def read(name):
    with open(pm + name) as f:
        txt = f.read()
    if re.search(r'format\s+binary', txt[:2000]):
        sys.exit(f'{name} is binary: set writeFormat ascii')
    i = txt.index('\n(', txt.index('FoamFile'))
    return txt[:i + 2], txt[i + 2:]


m = re.search(r'\b' + patch + r'\s*\{([^}]*)\}', read('boundary')[1])
if not m:
    sys.exit(f'no patch {patch}')
nFaces = int(re.search(r'nFaces\s+(\d+)', m.group(1)).group(1))
start = int(re.search(r'startFace\s+(\d+)', m.group(1)).group(1))

faces = re.findall(r'\d+\s*\(([^)]*)\)', read('faces')[1])
on = set()
for f in faces[start:start + nFaces]:
    on.update(int(v) for v in f.split())

head, body = read('points')
lines = body.split('\n')
dz = 0.0
for i in on:
    x, y, z = lines[i + 1].strip('()').split()
    dz = max(dz, abs(float(z)))
    lines[i + 1] = f'({x} {y} 0)'
with open(pm + 'points', 'w') as f:
    f.write(head + '\n'.join(lines))
print(f'{len(on)} points of {patch} set to z = 0, largest move {dz:.3g} m')
