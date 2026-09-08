#!/usr/bin/env python3
"""
blockTools.py -- a thin structured-block layer over the gmsh built-in kernel.

gmsh gives you points, curves, surfaces, volumes.  A multiblock hex mesh needs
them SHARED: two neighbouring blocks must reference the identical curve object,
or the transfinite meshes on either side will not match and the volumes will
not conform.  Everything here exists to make sharing automatic -- every entity
is memoised on its defining data, so asking for the same edge twice returns the
same tag rather than a duplicate sitting on top of it.

That single property is what makes a 30-block decomposition tractable by hand.
"""

import numpy as np
import gmsh


class Deck:
    """Memoising builder for a conforming structured multiblock model."""

    def __init__(self, model_name='mesh'):
        gmsh.initialize()
        gmsh.option.setNumber('General.Terminal', 0)

        gmsh.model.add(model_name)
        self.g = gmsh.model.geo
        self._pt, self._cv, self._sf, self._vl = {}, {}, {}, {}
        self._cv_ends = {}          # curve tag -> (start point, end point)
        self._cv_dir = {}           # curve tag -> direction class label
        self.surf_tags = {}         # physical name -> [surface tags]

    # ------------------------------------------------------------- points --
    def P(self, x, u, v):
        """Point at axial x, cross-plane (u, v).  Quarter model: y = u, z = -v."""
        k = (round(x, 10), round(u, 10), round(v, 10))
        if k not in self._pt:
            self._pt[k] = self.g.addPoint(k[0], k[1], -k[2])
        return self._pt[k]

    # ------------------------------------------------------------- curves --
    def _register(self, key, tag, p, q, nc):
        self._cv[key] = tag
        self._cv_ends[tag] = (p, q)
        self._cv_dir[tag] = nc          # (n_cells, progression coefficient)
        return tag

    def line(self, p, q, nc):
        key = ('L', min(p, q), max(p, q))
        if key in self._cv:
            return self._cv[key]
        return self._register(key, self.g.addLine(p, q), p, q, nc)

    def arc(self, p, q, x_centre, nc):
        """Circle arc about the x-axis, in the plane x = x_centre."""
        key = ('A', min(p, q), max(p, q))
        if key in self._cv:
            return self._cv[key]
        c = self.P(x_centre, 0.0, 0.0)
        return self._register(key, self.g.addCircleArc(p, c, q), p, q, nc)

    def spline(self, p, q, sampler, nc, n_ctrl=81):
        """Spline through `sampler(t)`, t in [0, 1], with the real progression.

        Note carefully how gmsh meshes this.  Transfinite points on a curve are
        spaced by ARC LENGTH, not by parameter -- so control points cannot be
        used to dictate node positions, and a first attempt at that produced a
        wall meridian with a uniform distribution sitting opposite a strongly
        graded straight edge.  The blocks sheared and non-orthogonality reached
        89.5 degrees.

        The right construction is the plain one: sample densely so the CURVE is
        accurate, hand gmsh the real progression, and rely on arc length and x
        agreeing to within sqrt(1 + r'^2) - 1.  On this body r' <= 0.2, so the
        two distributions differ by at most ~2 % of a block length, which is a
        fraction of one cell.
        """
        key = ('S', min(p, q), max(p, q))
        if key in self._cv:
            return self._cv[key]
        mids = [self.P(*sampler(float(t))) for t in np.linspace(0.0, 1.0, n_ctrl)[1:-1]]
        tag = self.g.addSpline([p] + mids + [q])
        return self._register(key, tag, p, q, nc)

    # ----------------------------------------------------------- surfaces --
    def face(self, corners, curves, physical=None):
        """Transfinite quad through 4 ordered corners and their 4 curves."""
        key = tuple(sorted(curves))
        if key in self._sf:
            tag = self._sf[key]
        else:
            signed = []
            for c, a in zip(curves, corners):
                s, _ = self._cv_ends[c]
                signed.append(c if s == a else -c)
            loop = self.g.addCurveLoop(signed)
            tag = self.g.addSurfaceFilling([loop])
            self.g.mesh.setTransfiniteSurface(tag, 'Left', list(corners))
            self.g.mesh.setRecombine(2, tag)
            self._sf[key] = tag
        if physical:
            self.surf_tags.setdefault(physical, []).append(tag)
        return tag

    # ------------------------------------------------------------ volumes --
    def block(self, faces, corners8):
        key = tuple(sorted(faces))
        if key in self._vl:
            return self._vl[key]
        loop = self.g.addSurfaceLoop(list(faces))
        tag = self.g.addVolume([loop])
        self.g.mesh.setTransfiniteVolume(tag, list(corners8))
        self._vl[key] = tag
        return tag

    # --------------------------------------------------------- transfinite --
    def apply_transfinite(self):
        for tag, (n, c) in self._cv_dir.items():
            self.g.mesh.setTransfiniteCurve(tag, int(n) + 1, 'Progression', float(c))

    def stats(self):
        return dict(points=len(self._pt), curves=len(self._cv),
                    surfaces=len(self._sf), volumes=len(self._vl))
