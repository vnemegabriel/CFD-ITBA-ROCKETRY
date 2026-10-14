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
(four, at finRoll + 0/90/180/270 deg, biconvex circular-arc section of constant
thickness, LE, TE and tip rounded with finEdgeRadius).  The shoulder and the
base edge are rounded with shoulderRadius and baseEdgeRadius, and a fillet of
finRootRadius runs all round each fin root.  Nose and body+fins+boattail are
two closed shells that meet on coplanar caps with rings offset half a facet.
All lengths in metres.
Nothing is written unless the parameters are consistent and every shell is
closed and outward-facing.

Next to out.stl go out_finEdges.obj, the leading and trailing edges of the
fins as lines, and out_finRoots.obj, the fin roots, for edgeMeshRefinement in
meshDict.
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
p.add_argument('--finRoll',       type=float, default=0.0,     help='angle of the first fin about the axis, deg')
p.add_argument('--finRootRadius', type=float, default=0.003,   help='fillet where the fin meets the body, all round the root')
p.add_argument('--shoulderRadius', type=float, default=0.003,  help='round between the cylinder and the boattail')
p.add_argument('--baseEdgeRadius', type=float, default=0.003,  help='round between the boattail and the base')
p.add_argument('--nFillet',       type=int,   default=6,       help='facets on a quarter circle of the fillet and the rounds')

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
rf, rs, rb = a.finRootRadius, a.shoulderRadius, a.baseEdgeRadius
thTail = math.atan2(R - Rb, a.tailLength)
L1 = rs * math.tan(thTail / 2)
L2 = rb * math.tan((0.5 * math.pi - thTail) / 2)

errors = [m for bad, m in [
    (min(rf, rs, rb) < 0,            'finRootRadius, shoulderRadius and baseEdgeRadius must be >= 0'),
    (rf >= a.finSpan / a.nSpan,      'finRootRadius must be below the first span station'),
    (rf >= xRootLE - Ln,             'finRootRadius runs onto the nose'),
    (L1 + L2 >= a.tailLength / math.cos(thTail), 'shoulderRadius and baseEdgeRadius do not fit on the boattail'),
    (L2 >= Rb,                       'baseEdgeRadius too large for the base'),
    (xRootTE + rf >= xBase - L2 * math.cos(thTail), 'finRootRadius runs onto the base edge'),
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


roll = math.radians(a.finRoll)


def rBody(x):
    """Radius of the body of revolution from the nose to the start of the
    base edge round: cylinder, shoulder round, cone."""
    if x <= xTail - L1:
        return R
    if x <= xTail + L1 * math.cos(thTail):
        return R - rs + math.sqrt(max(rs * rs - (x - xTail + L1) ** 2, 0.0))
    return R - (x - xTail) * math.tan(thTail)


def arc(cx, cr, rad, q0, q1):
    """Points (x, r) on a circle, q0 to q1 inclusive, about nFillet per quarter."""
    n = max(1, math.ceil(a.nFillet * abs(q1 - q0) / (0.5 * math.pi)))
    return [(cx + rad * math.cos(q0 + (q1 - q0) * j / n), cr + rad * math.sin(q0 + (q1 - q0) * j / n))
            for j in range(n + 1)]


# The profile aft of the shoulder: shoulder round, cone, base edge round, base.
aftProfile = []
if rs > 0:
    aftProfile += arc(xTail - L1, R - rs, rs, 0.5 * math.pi, 0.5 * math.pi - thTail)
else:
    aftProfile.append((xTail, R))
if rb > 0:
    q0 = math.atan2(Rb + L2 * math.sin(thTail) - (Rb - L2), xBase - L2 * math.cos(thTail) - (xBase - rb))
    aftProfile += arc(xBase - rb, Rb - L2, rb, q0, 0.0)
else:
    aftProfile.append((xBase, Rb))


def ring(s, x, r, offset):
    offset += roll / (2 * math.pi) * a.nTheta
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

# Body, boattail and fins, one shell.  Each fin stands on a fillet: a ring of
# arcs from the fin surface at height finRootRadius down to the body, a
# finRootRadius off the root contour.  The body between two fins is a sector
# of nq columns bounded by the outer edges of the two fillets.
body = Solid('body')
nq = a.nTheta // 4
root = section(a.finRootChord)
nc = len(root) // 2
xFoot = [xRootLE + q[0] + rf * q[2] for q in root[:nc + 1]]
zFoot = [q[1] + rf * q[3] for q in root[:nc + 1]]
fF = rf / a.finSpan
secF = section(a.finRootChord + fF * (a.finTipChord - a.finRootChord))
fins = []
for k in range(4):
    phi = roll + k * 0.5 * math.pi
    er, et = (math.cos(phi), math.sin(phi)), (-math.sin(phi), math.cos(phi))

    def at(x, r, z):
        return body.pt(x, r * er[0] + z * et[0], r * er[1] + z * et[1])

    # the fillet, from the body (q = pi/2) up to the fin (q = 0), in (x, r, z)
    # about this fin: P = K + (F - K)(1 - sin q) + (B - K)(1 - cos q), K on the
    # root contour, F on the fin at height rf, B on the body
    stations = []
    if rf > 0:
        arcs = []
        for (x, z, nx, nz), (xf, zf, _, _) in zip(root, secF):
            K = (xRootLE + x, math.sqrt(rBody(xRootLE + x) ** 2 - z * z), z)
            F = (xRootLE + fF * a.finSweep + xf, R + rf, zf)
            xb, zb = xRootLE + x + rf * nx, z + rf * nz
            B = (xb, math.sqrt(rBody(xb) ** 2 - zb * zb), zb)
            arcs.append((K, F, B))
        for n in range(a.nFillet, -1, -1):
            q = 0.5 * math.pi * n / a.nFillet
            stations.append([at(*(K[d] + (F[d] - K[d]) * (1 - math.sin(q)) + (B[d] - K[d]) * (1 - math.cos(q))
                                  for d in range(3))) for K, F, B in arcs])
    rho = a.finEdgeRadius
    for i in range(1 if rf > 0 else 0, a.nSpan + a.nRound):
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
            pts.append(at(xLE + x, math.sqrt(R * R - z * z) if i == 0 else r, z))
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
    x, r = xFoot[i], rBody(xFoot[i])
    w = math.asin(zFoot[i] / r)
    col = []
    for k in range(4):
        phi0 = roll + k * 0.5 * math.pi + w
        dphi = (0.5 * math.pi - 2 * w) / nq
        col.append([upper(k, i)] + [body.pt(x, r * math.cos(phi0 + j * dphi), r * math.sin(phi0 + j * dphi))
                                    for j in range(1, nq)] + [lower((k + 1) % 4, i)])
    columns.append(col)
ringLE = [v for col in columns[0] for v in col[:-1]]
ringTE = [v for col in columns[-1] for v in col[:-1]]
r0 = ring(body, xBody, R, 0.0)
strip(body, r0, ringLE, 'body')
for i in range(nc):
    name = 'boattail' if xFoot[i] >= xTail - 1e-9 else 'body'
    for k in range(4):
        for j in range(nq):
            body.quad(columns[i][k][j], columns[i + 1][k][j], columns[i + 1][k][j + 1], columns[i][k][j + 1], name)
last = ringTE
for x, r in aftProfile:
    if x > xFoot[-1] + 1e-9:
        rg = ring(body, x, r, 0.0)
        strip(body, last, rg, 'boattail' if x > xTail else 'body')
        last = rg
cap(body, r0, False, 'body')
cap(body, last, True, 'boattail')

solids = [nose, body]

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

def writeLines(suffix, segments):
    """segments: ((x0, r0), (x1, r1)) in the plane of each fin."""
    path = (a.out[:-4] if a.out.endswith('.stl') else a.out) + suffix
    with open(path, 'w') as out:
        n = 0
        for k in range(4):
            phi = roll + k * 0.5 * math.pi
            c, s_ = math.cos(phi), math.sin(phi)
            for (x0, r0), (x1, r1) in segments:
                out.write('v %.9g %.9g %.9g\nv %.9g %.9g %.9g\n' % (x0, r0 * c, r0 * s_, x1, r1 * c, r1 * s_))
                out.write(f'l {n + 1} {n + 2}\n')
                n += 2
    return path


edges = writeLines('_finEdges.obj', (((xRootLE, R), (xRootLE + a.finSweep, rTip)),
                                     ((xRootTE, R), (xRootLE + a.finSweep + a.finTipChord, rTip))))
roots = writeLines('_finRoots.obj', (((xRootLE, R), (xRootTE, R)),))

print(f'{a.out} + {edges} + {roots}: ' + ', '.join(f'{n} {len(t)} tris' for n, t in groups.items()))
print(f'length {xBase:.5f}  fins x {xRootLE:.5f}..{xRootLE + a.finSweep + a.finTipChord:.5f}  '
      f'r {R:.5f}..{rTip:.5f}  t {t * 1e3:.2f} mm  edge radius {a.finEdgeRadius * 1e3:.2f} mm')
