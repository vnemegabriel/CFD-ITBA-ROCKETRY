#!/usr/bin/env python3
"""
buildHexBody.py -- structured all-hexahedral mesh of the Aconcagua body.

TOPOLOGY, AND WHY IT IS THIS ONE
================================
The fluid around a pointed body of revolution is topologically an annulus
whose hole CLOSES at the nose tip and REOPENS at the base.  A single radial
O-grid cannot cover it: carried into the tip, every radial line converges on
one point and the cells collapse.  That is a property of the shape, not of the
mesher, and it is why a naive hex grid on a rocket always fails in exactly two
places.

So the topology changes precisely where the geometry forces it to, and nowhere
else.  The quarter cross-plane is tiled by five levels:

    level 0  core   -- butterfly square       } present only where the
    level 1  ring   -- square -> inner circle }   AXIS IS FLUID
    level 2  shell 1  boundary layer and near wake
    level 3  shell 2  mid field
    level 4  shell 3  out to the farfield

Upstream of the butterfly cap and downstream of the base the axis is fluid, so
core and ring exist; between them the axis is inside the solid and only the
annuli do.  The upstream core/ring blocks terminate ON the nose cap; the
downstream ones begin ON the flat base.  No cell is degenerate, and there is
no singular line anywhere in the fluid.

WHAT THIS BUYS OVER snappyHexMesh
=================================
The wall is a BOUNDARY of the block structure, not a TARGET for a projection.
A wall node is placed at r_body(x) by evaluating a formula; it does not migrate
there over nSolveIter iterations under a displacement scaling that can -- and
on this case did -- get "set to 0".  The boundary layer is the first ~31 cells
of a radial index, not the residue of a shrink-and-grow that decayed from
75.9 % coverage to 0.47 %.  There is nothing to snap, so nothing to control.
"""

import os
import numpy as np
import gmsh

import aconcaguaGeom as G
import meshParams as MP
from blockTools import Deck

SQ2 = np.sqrt(2.0)


def core_edge(M, k):
    """The core boundary edge facing azimuthal block k, as (start, end).

    Walking anticlockwise from theta = 0: up the u = a side for k < M, then back
    along v = a to the axis.  A SINGLE quad core can only present two edges,
    which is exactly what caps a butterfly at two azimuthal blocks per quadrant
    -- hence the M x M grid.
    """
    if k < M:
        return ('c', M, k), ('c', M, k + 1)
    m = k - M
    return ('c', M - m, M), ('c', M - m - 1, M)


def cross_quads(nlev, has_core, N):
    """The cross-plane block tiling, as node keys.  One truth.

        ('c', i, j)    butterfly core grid node, i, j in [0, M],  M = N/2
        ('r', lev, j)  ring / annulus node, j in [0, N] around the quadrant
    """
    M = N // 2
    q = []
    if has_core:
        for i in range(M):
            for j in range(M):
                q.append((('core', i, j),
                          [('c', i, j), ('c', i + 1, j),
                           ('c', i + 1, j + 1), ('c', i, j + 1)]))
        for k in range(N):
            a, b = core_edge(M, k)
            q.append((('ring', k), [a, ('r', 1, k), ('r', 1, k + 1), b]))
    for lev in range(1, nlev - 1):
        for k in range(N):
            q.append((('ann', lev, k),
                      [('r', lev, k), ('r', lev + 1, k),
                       ('r', lev + 1, k + 1), ('r', lev, k + 1)]))
    return q


def cap_x(u, v):
    """x of the nose point directly under cross-plane position (u, v).

    The butterfly is pushed onto the nose ALONG X: (u, v) is preserved and only
    x moves.  r_nose is monotone, so the map is a bijection.
    """
    return float(G.x_nose_of_r(np.hypot(u, v)))


class Station:
    """One cross-plane cut of the block structure."""

    def __init__(self, name, x, ri, has_core, kind, shells, shell1_c, th):
        self.name, self.x, self.ri = name, float(x), float(ri)
        self.has_core, self.kind, self.shell1_c = has_core, kind, shell1_c
        self.th = th
        self.N = len(th) - 1
        self.M = self.N // 2
        self.a = MP.CORE_FRAC * self.ri
        # one zone in the stack opens out downstream of the base so the fine
        # radial band tracks the spreading wake -- MP.zone_r_out knows which
        self.zr = [MP.zone_r_out(k, self.x) for k in range(len(shells))]
        self.r1_out = self.zr[0]
        self.radii = [None, self.ri] + self.zr
        if has_core:
            assert self.a * SQ2 < self.ri, f'{name}: core square exceeds inner circle'

    def uv(self, key):
        if key[0] == 'c':
            _, i, j = key
            return (self.a * i / self.M, self.a * j / self.M)
        _, lev, j = key
        r = self.radii[lev]
        return (r * np.cos(self.th[j]), r * np.sin(self.th[j]))

    def on_wall(self, key):
        return self.kind == 'cap' and (key[0] == 'c' or key[1] <= 1)

    def pos(self, key):
        u, v = self.uv(key)
        return (cap_x(u, v) if self.on_wall(key) else self.x, u, v)


class BodyMesh:

    def __init__(self):
        self.d = MP.derived()
        self.th = MP.az_angles()
        self.N = MP.N_AZ_BLOCKS
        self.M = self.N // 2
        self.nlev = len(self.d['shells']) + 2
        self.D = Deck('aconcagua')
        self.wall_lateral, self.wall_cap = set(), set()
        self._layout()
        self.cnt = self._counts()

    # ================================================== segments & stations ==
    def _layout(self):
        """Expand the seven named segments into the actual block columns."""
        d = self.d
        ax = MP.axial_plan(d)
        segs = []
        for nm, _, _ in MP.segment_bounds(d):
            if nm == 'up':
                # the upstream blend surface is curved, so subdivide it for the
                # same reason as the nose
                for i, p in enumerate(MP.subdivide(
                        ax['up'], MP.up_stations(d, ax['up']))):
                    segs.append(dict(name=f'up{i}', **p, core=True,
                                     wall=None, inner='up'))
                segs[-1]['name'] = 'up_last'
            elif nm == 'nose':
                for i, p in enumerate(MP.subdivide(ax['nose'],
                                                   MP.nose_stations(d, ax['nose']))):
                    segs.append(dict(name=f'nose{i}', **p, core=False,
                                     wall=MP.PATCHES['nose'], inner='body'))
            elif nm.startswith('wake'):
                segs.append(dict(name=nm, **{k: ax[nm][k] for k in
                                             ('x0', 'x1', 'n', 'c')},
                                 core=True, wall=None, inner='wake'))
            else:                                  # cyl, cylfin, tail
                wall = MP.PATCHES['tail'] if nm == 'tail' else MP.PATCHES['body']
                segs.append(dict(name=nm, **{k: ax[nm][k] for k in
                                             ('x0', 'x1', 'n', 'c')},
                                 core=False, wall=wall, inner='body'))
        self.segs = segs

        # --- one station per block interface -------------------------------
        self.S, self.order = {}, []
        c_wall = d['shells'][0]['c']
        xs = [segs[0]['x0']] + [s['x1'] for s in segs]
        for i, x in enumerate(xs):
            nm = f'st{i}'
            kind = 'cap' if abs(x - d['x_cap']) < 1e-12 else 'plane'
            ri, core, bias = self._station_spec(x, kind, c_wall)
            self.S[nm] = Station(nm, x, ri, core, kind, d['shells'], bias,
                                 self.th)
            self.order.append(nm)
        for i, s in enumerate(segs):
            s['s0'], s['s1'] = self.order[i], self.order[i + 1]

    def _station_spec(self, x, kind, c_wall):
        """Inner radius, whether the axis is fluid, and the shell-1 bias."""
        d = self.d
        if x < d['x_cap'] - 1e-12:                       # upstream of the cap
            phi = (x - d['x_in']) / (d['x_cap'] - d['x_in'])
            ri = self.inner_radius('up', x)
            return ri, True, 1.0 + phi ** 3 * (c_wall - 1.0)
        if x <= G.X_BASE + 1e-12:                        # on the body
            return float(G.r_body(x)), abs(x - G.X_BASE) < 1e-12, c_wall
        phi = (x - G.X_BASE) / (d['x_out'] - G.X_BASE)   # in the wake
        return self.inner_radius('wake', x), True, 1.0 + (1.0 - phi) ** 2 * (c_wall - 1.0)

    def inner_radius(self, kind, x):
        """Radius of shell 1's inner boundary.  On the body this IS the wall."""
        d = self.d
        if kind == 'body':
            return float(G.r_body(x))
        if kind == 'up':
            return float(np.sqrt(np.maximum(self._blend(d['x_cap'] - x), 1e-24)))
        xi = (x - G.X_BASE) / (d['x_out'] - G.X_BASE)
        return G.R_BASE + (MP.r_wake_out() - G.R_BASE) * xi

    def _blend(self, t):
        """Ri^2 along the upstream blend; tangent to the wall at the cap rim."""
        d = self.d
        rc, sc_ = d['r_cap'], float(G.drdx_body(d['x_cap']))
        D = d['x_cap'] - d['x_in']
        lam = (MP.r_inlet() ** 2 - rc ** 2 - 2.0 * rc * sc_ * D) / D ** 2
        return rc ** 2 + 2.0 * rc * sc_ * t + lam * t * t

    # ============================================================== counts ===
    def _counts(self):
        coefs = MP.az_coefs()
        nc = MP.n_az_cells()
        c = {'ring': (MP.n_ring(), 1.0)}
        for k, co in enumerate(coefs):
            c[f'azv{k}'] = (nc, co)                    # arcs, and core v-edges
        for i in range(self.M):
            c[f'azu{i}'] = (nc, 1.0 / coefs[self.N - 1 - i])   # core u-edges
        for seg in self.segs:
            c[f"ax_{seg['name']}"] = (seg['n'], seg['c'])
        for k, sh in enumerate(self.d['shells'], start=1):
            for nm, st in self.S.items():
                if k == 1:
                    span = st.zr[0] - st.ri
                    h1 = ((sh['r_out'] - G.R_BODY) * (sh['c'] - 1.0)
                          / (sh['c'] ** sh['n'] - 1.0))
                    coef = MP.refit(span, sh['n'], min(h1, 0.9 * span / sh['n'])) \
                        if st.shell1_c != 1.0 else 1.0
                    c[f'shell{k}@{nm}'] = (sh['n'], coef)
                else:
                    c[f'shell{k}@{nm}'] = (sh['n'], sh['c'])
        return c

    # ============================================================== curves ===
    def _ray(self, st, ka, kb):
        u0, v0 = st.uv(ka)
        u1, v1 = st.uv(kb)
        def f(t):
            u, v = u0 + t * (u1 - u0), v0 + t * (v1 - v0)
            return (cap_x(u, v), u, v)
        return f

    @staticmethod
    def _az(ka, kb):
        """Canonical direction and family for an azimuthal curve.

        Curves are always CREATED in the canonical direction (increasing index)
        even when a block asks for them backwards, because a progression is
        meaningless without one; face() recovers the sign from the endpoints.

        `azv{j}` is the arc of block j AND the core's v-direction edges -- the
        same edge of the same ring block, so they must share a distribution.
        `azu{i}` is the core's u-direction edge, which faces block N-1-i and is
        traversed the other way, hence the reciprocal coefficient in _counts.
        """
        if ka[0] == 'r':
            j = min(ka[2], kb[2])
            return ('r', ka[1], j), ('r', ka[1], j + 1), f'azv{j}'
        _, i0, j0 = ka
        _, i1, j1 = kb
        if i0 == i1:
            j = min(j0, j1)
            return ('c', i0, j), ('c', i0, j + 1), f'azv{j}'
        i = min(i0, i1)
        return ('c', i, j0), ('c', i + 1, j0), f'azu{i}'

    def c_cross(self, st, ka, kb):
        D = self.D
        same = (ka[0] == 'c' and kb[0] == 'c') or \
               (ka[0] == 'r' and kb[0] == 'r' and ka[1] == kb[1])
        if same:                                             # azimuthal family
            ka, kb, fam = self._az(ka, kb)
            p, q = D.P(*st.pos(ka)), D.P(*st.pos(kb))
            if ka[0] == 'r':
                return D.arc(p, q, st.x, self.cnt[fam])      # an exact circle
            if st.on_wall(ka):
                return D.spline(p, q, self._ray(st, ka, kb), self.cnt[fam])
            return D.line(p, q, self.cnt[fam])
        p, q = D.P(*st.pos(ka)), D.P(*st.pos(kb))
        # RADIAL FAMILY.  Canonicalise inward -> outward before creating the
        # curve.  A block asks for this edge forwards on one side and BACKWARDS
        # on the other (it is edge 0 of one annulus block and edge 2 of its
        # neighbour), and whichever call lands first fixes the direction the
        # progression runs in.  Left uncanonicalised, column 0 got the intended
        # wall-clustered stack and every other column got its reciprocal --
        # first radial cell 11.5 mm against 30.5 mm at the same station.  That
        # was the real source of the azimuthal non-orthogonality, not the
        # arc-to-chord error it was mistaken for.
        lo, hi = sorted([ka, kb],
                        key=lambda k: (0, 0) if k[0] == 'c' else (1, k[1]))
        lab = 'ring' if lo[0] == 'c' else f'shell{hi[1] - 1}@{st.name}'
        p, q = D.P(*st.pos(lo)), D.P(*st.pos(hi))
        if st.on_wall(hi):
            return D.spline(p, q, self._ray(st, lo, hi), self.cnt[lab])
        return D.line(p, q, self.cnt[lab])

    def c_axial(self, seg, s0, s1, key):
        D = self.D
        p, q = D.P(*s0.pos(key)), D.P(*s1.pos(key))
        lab = f"ax_{seg['name']}"
        straight = (seg['inner'] == 'body'
                    and seg['x0'] >= G.X_BODY_1 - 1e-12)   # cylinder + boattail
        if key[0] == 'r' and key[1] == 1 and not straight \
                and seg['inner'] in ('body', 'up'):
            th = float(self.th[key[2]])
            kind, x0, dx = seg['inner'], s0.x, s1.x - s0.x
            def f(t):
                x = x0 + t * dx
                r = self.inner_radius(kind, x)
                return (x, r * np.cos(th), r * np.sin(th))
            return D.spline(p, q, f, self.cnt[lab])
        return D.line(p, q, self.cnt[lab])

    # =============================================================== faces ===
    def face_cross(self, st, quad, physical=None):
        cor = [self.D.P(*st.pos(k)) for k in quad]
        cvs = [self.c_cross(st, quad[i], quad[(i + 1) % 4]) for i in range(4)]
        return self.D.face(cor, cvs, physical)

    def face_side(self, seg, s0, s1, ka, kb, physical=None):
        D = self.D
        cor = [D.P(*s0.pos(ka)), D.P(*s0.pos(kb)),
               D.P(*s1.pos(kb)), D.P(*s1.pos(ka))]
        cvs = [self.c_cross(s0, ka, kb), self.c_axial(seg, s0, s1, kb),
               self.c_cross(s1, ka, kb), self.c_axial(seg, s0, s1, ka)]
        return D.face(cor, cvs, physical)

    def _side_patch(self, seg, st, ka, kb):
        P = MP.PATCHES
        (ua, va), (ub, vb) = st.uv(ka), st.uv(kb)
        if abs(va) < 1e-12 and abs(vb) < 1e-12:
            return P['symmetry']
        if abs(ua) < 1e-12 and abs(ub) < 1e-12:
            return P['symmetry']
        rr = ka[0] == 'r' and kb[0] == 'r' and ka[1] == kb[1]
        if rr and ka[1] == self.nlev - 1:
            return P['farfield']                       # outermost cylinder
        if rr and ka[1] == 1 and seg['wall']:
            return seg['wall']                         # the body wall
        return None

    # =============================================================== build ===
    def build(self):
        P = MP.PATCHES
        first, last = self.order[0], self.order[-1]
        for seg in self.segs:
            s0, s1 = self.S[seg['s0']], self.S[seg['s1']]
            for name, quad in cross_quads(self.nlev, seg['core'], self.N):
                inner = name[0] in ('core', 'ring')
                ph0 = P['inlet'] if seg['s0'] == first else None
                ph1 = P['outlet'] if seg['s1'] == last else None
                if inner and s1.kind == 'cap':
                    ph1 = P['nose']                     # the nose cap itself
                    self.wall_cap.add(self.face_cross(s1, quad, None))
                if inner and abs(s0.x - G.X_BASE) < 1e-12:
                    ph0 = P['tail']                     # the flat base disc
                f0 = self.face_cross(s0, quad, ph0)
                f1 = self.face_cross(s1, quad, ph1)
                sides = []
                for i in range(4):
                    ka, kb = quad[i], quad[(i + 1) % 4]
                    sf = self.face_side(seg, s0, s1, ka, kb,
                                        self._side_patch(seg, s0, ka, kb))
                    if ka[0] == 'r' and kb[0] == 'r' and ka[1] == kb[1] == 1 \
                            and seg['wall']:
                        self.wall_lateral.add(sf)
                    sides.append(sf)
                corners = ([self.D.P(*s0.pos(k)) for k in quad] +
                           [self.D.P(*s1.pos(k)) for k in quad])
                self.D.block([f0, f1] + sides, corners)
        self.D.apply_transfinite()
        self.D.g.synchronize()
        for nm, tags in self.D.surf_tags.items():
            g = gmsh.model.addPhysicalGroup(2, sorted(set(tags)))
            gmsh.model.setPhysicalName(2, g, nm)
        g = gmsh.model.addPhysicalGroup(3, [v for _, v in gmsh.model.getEntities(3)])
        gmsh.model.setPhysicalName(3, g, 'internal')


if __name__ == '__main__':
    m = BodyMesh()
    m.build()
    print('geometry:', m.D.stats(), ' blocks:', len(m.segs))
    gmsh.finalize()
