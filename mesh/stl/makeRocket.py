#!/usr/bin/env python3
"""Writes a parametric Aconcagua-type rocket as an ASCII STL.

    ./makeRocket.py out.stl                       the current Aconcagua
    ./makeRocket.py long.stl --bodyLength 2.5
    ./makeRocket.py big.stl --diameter 0.2 --finSpan 0.2 --finFromBase 0.2
    ./makeRocket.py new.stl --dims new.dims       dimensions from a file
    ./makeRocket.py new.stl --dims new.dims --finSpan 0.18

A dims file holds one `name value` per line (a trailing `;` and `#` or `//`
comments are allowed); names are the options below without `--`.  Options on
the command line win over the file, the file over the defaults.

Solids: nosecone (von Karman), body (cylinder), boattail (cone + base), fins
(four, at 0/90/180/270 deg).  All lengths in metres.  Nothing is written unless
the parameters are consistent and the surface is closed and outward-facing.
"""
import argparse, math, re, sys
from collections import Counter

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('out')
p.add_argument('--dims', help='file of `name value` lines')
p.add_argument('--diameter',      type=float, default=0.151)
p.add_argument('--noseLD',        type=float, default=0.8 / 0.151, help='nose length / diameter')
p.add_argument('--bodyLength',    type=float, default=2.03,  help='cylinder only')
p.add_argument('--tailLength',    type=float, default=0.125)
p.add_argument('--baseDiameter',  type=float, default=0.110)
p.add_argument('--finRootChord',  type=float, default=0.30047)
p.add_argument('--finTipChord',   type=float, default=0.15)
p.add_argument('--finSpan',       type=float, default=0.16,    help='from the body surface to the tip')
p.add_argument('--finSweep',      type=float, default=0.25078, help='tip LE aft of root LE')
p.add_argument('--finFromBase',   type=float, default=0.12531, help='root TE to base plane')
p.add_argument('--finThickness',  type=float, default=0.006)
p.add_argument('--finLEBevel',    type=float, default=0.03316)
p.add_argument('--finTEBevel',    type=float, default=0.02118)
p.add_argument('--finEmbed',      type=float, default=0.0005,  help='root sunk into the body')
p.add_argument('--nTheta',        type=int,   default=96)
p.add_argument('--nNose',         type=int,   default=60)

def read_dims(path):
    types = {x.dest: x.type for x in p._actions if x.type in (float, int)}
    vals = {}
    for n, line in enumerate(open(path), 1):
        line = re.split(r'#|//', line)[0].strip().rstrip(';').strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) != 2 or parts[0] not in types:
            sys.exit(f'{path}:{n}: expected `name value` with name in {", ".join(types)}')
        try:
            vals[parts[0]] = types[parts[0]](parts[1])
        except ValueError:
            sys.exit(f'{path}:{n}: {parts[0]} needs a number, got {parts[1]}')
    return vals

pre, _ = p.parse_known_args()
if pre.dims:
    p.set_defaults(**read_dims(pre.dims))
a = p.parse_args()

R, Rb = a.diameter / 2, a.baseDiameter / 2
Ln = a.noseLD * a.diameter
x1, x2 = Ln, Ln + a.bodyLength
xb = x2 + a.tailLength
t2 = a.finThickness / 2
rootR = R - a.finEmbed
tipY = R + a.finSpan
rootTE = xb - a.finFromBase
rootLE = rootTE - a.finRootChord
tipLE = rootLE + a.finSweep
tipTE = tipLE + a.finTipChord

errors = []
def need(ok, msg):
    if not ok: errors.append(msg)

for k in ('diameter', 'noseLD', 'bodyLength', 'tailLength', 'baseDiameter', 'finRootChord',
          'finTipChord', 'finSpan', 'finThickness', 'finLEBevel', 'finTEBevel'):
    need(getattr(a, k) > 0, f'--{k} must be > 0')
need(Rb <= R, '--baseDiameter must not exceed --diameter')
need(a.nTheta >= 16 and a.nTheta % 4 == 0, '--nTheta must be a multiple of 4, >= 16 (fins sit on grid angles)')
need(a.nNose >= 10, '--nNose must be >= 10')
need(0 < a.finEmbed < R / 4, '--finEmbed must be > 0 (the fin must cut into the body) and small')
need(rootLE >= x1, f'fin root LE x={rootLE:.4f} is on the nose (starts at {x1:.4f}): '
                   'lower --finFromBase or shorten --finRootChord')
need(rootTE <= x2, f'fin root TE x={rootTE:.4f} is on the boattail (starts at {x2:.4f}): '
                   'the root would float off the cone; raise --finFromBase')
need(a.finLEBevel + a.finTEBevel < min(a.finRootChord, a.finTipChord),
     'LE + TE bevels must be shorter than both chords')
need(t2 < rootR * math.sin(math.pi / 8), 'fins too thick for four around this body')
if errors:
    sys.exit('refused:\n  ' + '\n  '.join(errors))

tris = {'nosecone': [], 'body': [], 'boattail': [], 'fins': []}
N = a.nTheta
ang = [2 * math.pi * k / N for k in range(N)]

def ring(x, r):
    return [(x, r * math.cos(t), r * math.sin(t)) for t in ang]

def band(solid, A, B):
    for k in range(N):
        k1 = (k + 1) % N
        tris[solid] += [(A[k], B[k1], B[k]), (A[k], A[k1], B[k1])]

def r_nose(x):
    th = math.acos(max(-1.0, min(1.0, 1 - 2 * x / Ln)))
    return R / math.sqrt(math.pi) * math.sqrt(max(th - math.sin(2 * th) / 2, 0.0))

xs = [Ln * (i / a.nNose) ** 2 for i in range(1, a.nNose + 1)]
rings = [ring(x, r_nose(x)) for x in xs]
rings[-1] = ring(x1, R)
tip = (0.0, 0.0, 0.0)
for k in range(N):
    tris['nosecone'].append((tip, rings[0][(k + 1) % N], rings[0][k]))
for A, B in zip(rings, rings[1:]):
    band('nosecone', A, B)

band('body', ring(x1, R), ring(x2, R))
band('boattail', ring(x2, R), ring(xb, Rb))
base, centre = ring(xb, Rb), (xb, 0.0, 0.0)
for k in range(N):
    tris['boattail'].append((base[k], base[(k + 1) % N], centre))

def fin_section(xLE, xTE, y_of):
    xs = [xLE, xLE + a.finLEBevel, xTE - a.finTEBevel, xTE]
    ts = [0.0, t2, t2, 0.0]
    upper = [(x, y_of(t), t) for x, t in zip(xs, ts)]
    lower = [(x, y_of(t), -t) for x, t in zip(xs, ts)]
    return upper, lower

root_u, root_l = fin_section(rootLE, rootTE, lambda t: math.sqrt(rootR ** 2 - t ** 2))
tip_u, tip_l = fin_section(tipLE, tipTE, lambda t: tipY)
one = []
for i in range(3):
    one += [(root_u[i], root_u[i + 1], tip_u[i + 1]), (root_u[i], tip_u[i + 1], tip_u[i])]
    one += [(root_l[i], tip_l[i + 1], root_l[i + 1]), (root_l[i], tip_l[i], tip_l[i + 1])]
hexagon = lambda u, l: [u[0], u[1], u[2], u[3], l[2], l[1]]
rh, th_ = hexagon(root_u, root_l), hexagon(tip_u, tip_l)
for i in range(1, 5):
    one.append((rh[0], rh[i + 1], rh[i]))
    one.append((th_[0], th_[i], th_[i + 1]))
for q in range(4):
    c, s = math.cos(q * math.pi / 2), math.sin(q * math.pi / 2)
    rot = lambda v: (v[0], c * v[1] - s * v[2], s * v[1] + c * v[2])
    tris['fins'] += [tuple(rot(v) for v in t) for t in one]

def sub(u, v): return (u[0] - v[0], u[1] - v[1], u[2] - v[2])
def cross(u, v): return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
def dot(u, v): return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]

def check_closed(name, ts):
    key = lambda v: tuple(round(c, 12) for c in v)
    edges = Counter()
    vol = 0.0
    for t in ts:
        n = cross(sub(t[1], t[0]), sub(t[2], t[0]))
        if math.sqrt(dot(n, n)) < 1e-14:
            errors.append(f'{name}: degenerate triangle at {t[0]}')
        vol += dot(t[0], cross(t[1], t[2])) / 6
        for i in range(3):
            e = (key(t[i]), key(t[(i + 1) % 3]))
            edges[e] += 1
    open_edges = sum(1 for (u, v), n in edges.items() if n != 1 or edges[(v, u)] != 1)
    if open_edges:
        errors.append(f'{name}: {open_edges} edges not shared by exactly two consistently oriented triangles')
    if vol <= 0:
        errors.append(f'{name}: normals point inward (volume {vol:.3e})')
    return vol

vol_body = check_closed('body of revolution', tris['nosecone'] + tris['body'] + tris['boattail'])
for q in range(4):
    check_closed(f'fin {q}', tris['fins'][q * len(one):(q + 1) * len(one)])
if errors:
    sys.exit('refused, bad surface:\n  ' + '\n  '.join(errors[:10]))

with open(a.out, 'w') as f:
    for name, ts in tris.items():
        f.write(f'solid {name}\n')
        for t in ts:
            n = cross(sub(t[1], t[0]), sub(t[2], t[0]))
            m = math.sqrt(dot(n, n))
            f.write(' facet normal %.9g %.9g %.9g\n  outer loop\n' % tuple(c / m for c in n))
            for v in t:
                f.write('   vertex %.9g %.9g %.9g\n' % v)
            f.write('  endloop\n endfacet\n')
        f.write(f'endsolid {name}\n')

print(f'wrote {a.out}: {sum(map(len, tris.values()))} triangles, '
      f'L = {xb:.4f} m, D = {a.diameter:.4f} m, L/D = {xb / a.diameter:.2f}')
print('  nose 0 -> %.4f   body -> %.4f   boattail -> %.4f' % (x1, x2, xb))
print('  fins: root %.4f -> %.4f   tip %.4f -> %.4f   tip r %.4f' % (rootLE, rootTE, tipLE, tipTE, tipY))
print('meshDict boxes and wake cones are in absolute x: move them if the base moved '
      f'(base now at x = {xb:.4f}, Aconcagua 2.955)')
