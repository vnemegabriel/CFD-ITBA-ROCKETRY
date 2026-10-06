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
(four, at 0/90/180/270 deg, biconvex circular-arc section of constant
thickness, LE, TE and tip rounded with finEdgeRadius).  Nose, body+fins and boattail are three closed shells that meet
on coplanar caps with rings offset half a facet.  Body and fins are one
shell: the fin root is cut into the cylinder, so the fin-body junction is a
real edge of the surface.  All lengths in metres.
Nothing is written unless the parameters are consistent and every shell is
closed and outward-facing.

Next to out.stl goes out_finEdges.obj: the leading and trailing edges of the
fins as lines, for edgeMeshRefinement in meshDict.
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
p.add_argument('--finFromBase',   type=float, default=0.125,   help='root TE to base plane')
p.add_argument('--finThickness',  type=float, default=0.0117,  help='max thickness, at mid-chord, same at every span station')
p.add_argument('--nTheta',        type=int,   default=44,      help='facets around the axis')
p.add_argument('--nNose',         type=int,   default=60,      help='stations along the nose')
p.add_argument('--finEdgeRadius', type=float, default=0.001,   help='round on the LE, TE and tip, in the streamwise section')
p.add_argument('--nChord',        type=int,   default=24,      help='facets per fin side on the biconvex arc')
p.add_argument('--nRound',        type=int,   default=6,       help='facets on each edge round')
p.add_argument('--nSpan',         type=int,   default=8,       help='fin span stations')

pre, _ = p.parse_known_args()
if pre.dims:
    known = {a.dest: a.type for a in p._actions if a.type}
    vals = {}
    for n, line in enumerate(open(pre.dims), 1):
        line = re.split(r'#|//', line)[0].strip().rstrip(';').strip()
        if not line:
            continue
        name, value = line.split()
        if name not in known:
            sys.exit(f'{pre.dims}:{n}: unknown name {name}')
        vals[name] = known[name](value)
    p.set_defaults(**vals)
a = p.parse_args()

R, Rb = a.diameter / 2, a.baseDiameter / 2
Ln = a.noseLD * a.diameter
xBody, xTail = Ln, Ln + a.bodyLength
xBase = xTail + a.tailLength
t = a.finThickness
rTip = R + a.finSpan
xRootTE = xBase - a.finFromBase
xRootLE = xRootTE - a.finRootChord

errors = [m for bad, m in [
    (Rb >= R,                        'baseDiameter must be smaller than diameter'),
    (a.finTipChord > a.finRootChord, 'finTipChord larger than finRootChord'),
    (a.finFromBase < a.tailLength,   'fin root runs onto the boattail: finFromBase < tailLength'),
    (xRootLE < xBody,                'fin root starts on the nose'),
    (t >= R,                         'finThickness too large for the body'),
    (a.nTheta % 4 != 0,              'nTheta must be a multiple of 4'),
    (not 0 < a.finEdgeRadius < t / 2, 'finEdgeRadius must be between 0 and finThickness/2'),
    (t >= 0.5 * a.finTipChord,       'finThickness too large for the tip chord'),
] if bad]
if errors:
    sys.exit('\n'.join(errors))


class Solid:
    def __init__(self, name):
        self.name, self.pts, self.tris, self.names = name, [], [], []

    def pt(self, x, y, z):
        self.pts.append((x, y, z))
        return len(self.pts) - 1

    def quad(self, a, b, c, d, name=None):
        for tri in ((a, b, c), (a, c, d)):
            if len(set(tri)) == 3:
                self.tris.append(tri)
                self.names.append(name or self.name)

    def volume(self):
        v = 0.0
        for i, j, k in self.tris:
            (x1, y1, z1), (x2, y2, z2), (x3, y3, z3) = self.pts[i], self.pts[j], self.pts[k]
            v += x1 * (y2 * z3 - z2 * y3) - y1 * (x2 * z3 - z2 * x3) + z1 * (x2 * y3 - y2 * x3)
        return v / 6

    def check(self):
        directed = Counter((tri[e], tri[(e + 1) % 3]) for tri in self.tris for e in range(3))
        if any(c > 1 for c in directed.values()):
            sys.exit(f'{self.name}: inconsistent orientation')
        if any((b, a) not in directed for a, b in directed):
            sys.exit(f'{self.name}: surface not closed')
        if self.volume() < 0:
            self.tris = [(i, k, j) for i, j, k in self.tris]
        if self.volume() <= 0:
            sys.exit(f'{self.name}: zero volume')


def ring(s, x, r, offset):
    return [s.pt(x, r * math.cos(2 * math.pi * (j + offset) / a.nTheta),
                 r * math.sin(2 * math.pi * (j + offset) / a.nTheta)) for j in range(a.nTheta)]


def strip(s, r0, r1, name):
    for j in range(len(r0)):
        k = (j + 1) % len(r0)
        s.quad(r0[j], r1[j], r1[k], r0[k], name)


def cap(s, ring_, last, name):
    for j in range(1, len(ring_) - 1):
        s.tris.append((ring_[0], ring_[j + 1], ring_[j]) if last else (ring_[0], ring_[j], ring_[j + 1]))
        s.names.append(name)


def section(c):
    """Closed contour, LE at 0: upper LE->TE, then lower TE->LE (open), as
    (x, z, nx, nz) with the outward normal.  The convex hull of the biconvex
    arcs and a circle of finEdgeRadius at each edge: round, tangent segment,
    arc, tangent segment, round."""
    rho = a.finEdgeRadius
    Rc = (c * c + t * t) / (4 * t)
    d = Rc - 0.5 * t
    Dx, Dz = rho - 0.5 * c, d
    beta = math.atan2(Dz, Dx) - math.acos((Rc - rho) / math.hypot(Dx, Dz))
    nb = (math.cos(beta), math.sin(beta))
    T1 = (rho + rho * nb[0], rho * nb[1])
    T2 = (0.5 * c + Rc * nb[0], -d + Rc * nb[1])
    side = []
    for j in range(a.nRound):
        q = math.pi + (beta - math.pi) * j / a.nRound
        side.append((rho + rho * math.cos(q), rho * math.sin(q), math.cos(q), math.sin(q)))
    for j in range(a.nRound):
        u = j / a.nRound
        side.append((T1[0] + u * (T2[0] - T1[0]), T1[1] + u * (T2[1] - T1[1]), nb[0], nb[1]))
    q0, q1 = beta, math.pi - beta
    for j in range(a.nChord):
        q = q0 + (q1 - q0) * j / a.nChord
        side.append((0.5 * c + Rc * math.cos(q), -d + Rc * math.sin(q), math.cos(q), math.sin(q)))
    aft = [(c - x, z, -nx, nz) for x, z, nx, nz in reversed(side[:2 * a.nRound])]
    up = side + [(c - T2[0], T2[1], -nb[0], nb[1])] + aft[:-1] + [(c, 0.0, 1.0, 0.0)]
    up[0] = (0.0, 0.0, -1.0, 0.0)
    return up + [(x, -z, nx, -nz) for x, z, nx, nz in reversed(up)][1:-1]


nose = Solid('nosecone')
th = [math.pi * i / a.nNose for i in range(1, a.nNose)]
rings = [[nose.pt(0.0, 0.0, 0.0)] * a.nTheta]
rings += [ring(nose, Ln * (1 - math.cos(q)) / 2, R / math.sqrt(math.pi) * math.sqrt(q - math.sin(2 * q) / 2), 0.5) for q in th]
rings.append(ring(nose, Ln, R, 0.5))
for r0, r1 in zip(rings, rings[1:]):
    strip(nose, r0, r1, 'nosecone')
cap(nose, rings[-1], True, 'nosecone')

tail = Solid('boattail')
rA, rB = ring(tail, xTail, R, 0.5), ring(tail, xBase, Rb, 0.5)
strip(tail, rA, rB, 'boattail')
cap(tail, rA, False, 'boattail')
cap(tail, rB, True, 'boattail')

# Body and fins, one shell.  Each fin root contour lies on the cylinder; the
# body between two fins is a sector of nq columns bounded by those contours.
body = Solid('body')
nq = a.nTheta // 4
root = section(a.finRootChord)
nc = len(root) // 2
zRoot = [q[1] for q in root[:nc + 1]]
fins = []
for k in range(4):
    phi = k * 0.5 * math.pi
    er, et = (math.cos(phi), math.sin(phi)), (-math.sin(phi), math.cos(phi))
    stations = []
    rho = a.finEdgeRadius
    for i in range(a.nSpan + a.nRound):
        if i <= a.nSpan:
            f = i / a.nSpan * (a.finSpan - rho) / a.finSpan
            r, off = R + f * a.finSpan, 0.0
            xLE = xRootLE + f * a.finSweep
            sec = root if i == 0 else section(a.finRootChord + f * (a.finTipChord - a.finRootChord))
        else:
            q = 0.45 * math.pi * (i - a.nSpan) / (a.nRound - 1)
            r, off = rTip - rho + rho * math.sin(q), rho * (1 - math.cos(q))
        pts = []
        for x, z, nx, nz in sec:
            x, z = x - off * nx, z - off * nz
            rr = math.sqrt(R * R - z * z) if i == 0 else r
            pts.append(body.pt(xLE + x, rr * er[0] + z * et[0], rr * er[1] + z * et[1]))
        stations.append(pts)
    fins.append(stations)
    m = len(stations[0])
    for s0, s1 in zip(stations, stations[1:]):
        for j in range(m):
            jj = (j + 1) % m
            body.quad(s0[j], s1[j], s1[jj], s0[jj], 'fins')
    tip = stations[-1]
    cen = body.pt(*(sum(body.pts[v][d] for v in tip) / m for d in range(3)))
    for j in range(m):
        body.tris.append((cen, tip[(j + 1) % m], tip[j]))
        body.names.append('fins')


def upper(k, i):
    return fins[k][0][i]


def lower(k, i):
    return fins[k][0][(2 * nc - i) % (2 * nc)]


columns = []
for i in range(nc + 1):
    x = xRootLE + root[i][0]
    w = math.asin(zRoot[i] / R)
    col = []
    for k in range(4):
        phi0 = k * 0.5 * math.pi + w
        dphi = (0.5 * math.pi - 2 * w) / nq
        col.append([upper(k, i)] + [body.pt(x, R * math.cos(phi0 + j * dphi), R * math.sin(phi0 + j * dphi))
                                    for j in range(1, nq)] + [lower((k + 1) % 4, i)])
    columns.append(col)
ringLE = [v for col in columns[0] for v in col[:-1]]
ringTE = [v for col in columns[-1] for v in col[:-1]]
r0 = ring(body, xBody, R, 0.0)
r1 = ringTE if abs(xRootTE - xTail) < 1e-9 else ring(body, xTail, R, 0.0)
strip(body, r0, ringLE, 'body')
for c0, c1 in zip(columns, columns[1:]):
    for k in range(4):
        for j in range(nq):
            body.quad(c0[k][j], c1[k][j], c1[k][j + 1], c0[k][j + 1], 'body')
if r1 is not ringTE:
    strip(body, ringTE, r1, 'body')
cap(body, r0, False, 'body')
cap(body, r1, True, 'body')

solids = [nose, body, tail]

for s in solids:
    s.check()

groups = {n: [] for n in ('nosecone', 'body', 'boattail', 'fins')}
for sol in solids:
    for tri, name in zip(sol.tris, sol.names):
        groups[name].append([sol.pts[i] for i in tri])

with open(a.out, 'w') as out:
    for name, tris in groups.items():
        out.write(f'solid {name}\n')
        for P, Q, S in tris:
            u = [Q[d] - P[d] for d in range(3)]
            v = [S[d] - P[d] for d in range(3)]
            nrm = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
            ln = math.sqrt(sum(q * q for q in nrm)) or 1.0
            out.write(' facet normal %.9g %.9g %.9g\n  outer loop\n' % tuple(q / ln for q in nrm))
            for q in (P, Q, S):
                out.write('   vertex %.9g %.9g %.9g\n' % q)
            out.write('  endloop\n endfacet\n')
        out.write(f'endsolid {name}\n')

edges = a.out[:-4] + '_finEdges.obj' if a.out.endswith('.stl') else a.out + '_finEdges.obj'
with open(edges, 'w') as out:
    n = 0
    for phi in (0, 0.5 * math.pi, math.pi, 1.5 * math.pi):
        c, s_ = math.cos(phi), math.sin(phi)
        for x0, x1 in ((xRootLE, xRootLE + a.finSweep), (xRootTE, xRootLE + a.finSweep + a.finTipChord)):
            out.write('v %.9g %.9g %.9g\nv %.9g %.9g %.9g\n' % (x0, R * c, R * s_, x1, rTip * c, rTip * s_))
            out.write(f'l {n + 1} {n + 2}\n')
            n += 2

print(f'{a.out} + {edges}: ' + ', '.join(f'{n} {len(t)} tris' for n, t in groups.items()))
print(f'length {xBase:.5f}  fins x {xRootLE:.5f}..{xRootLE + a.finSweep + a.finTipChord:.5f}  '
      f'r {R:.5f}..{rTip:.5f}  t {t * 1e3:.2f} mm  edge radius {a.finEdgeRadius * 1e3:.2f} mm')
