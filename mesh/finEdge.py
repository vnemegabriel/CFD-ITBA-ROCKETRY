#!/usr/bin/env python3
"""
finEdge.py -- put the fin leading and trailing edges ON the mesh.

THE PROBLEM
===========
The fin planform edges are swept.  In the (x, r) plane of a symmetry plane the
leading edge runs at dx/dr = 1.5626 and the trailing edge at 0.6249, while the
grid there is x-stations crossed with r-shells.  The edges therefore cut across
the grid in diagonal, and finPatch.split_symm(), which can only call a WHOLE
face wall or symmetry, quantises the outline to the cell.

That rule is not a sloppy choice and swapping it for a centroid test does not
help: at assembly a symmetry-plane face can be merged into the interior only
if all four of its nodes are still at z = 0, so "any node displaced -> wall" is
forced.  The consequence is a patch that is the true planform DILATED by one
cell -- never smaller, always larger:

    H_SCALE 1   +10.3 % of planform area      50 span rows, 81 chord cells
    H_SCALE 3   +18.9 %                       16 span rows, 27 chord cells
    H_SCALE 6   +21.1 %                        8 span rows, 14 chord cells

and, because the wedge section has t -> 0 at the leading edge, the extra ring
of faces carries almost no thickness: a serrated flange of near-zero thickness
standing ahead of the real leading edge, which is the last thing a sharp
supersonic edge wants.

THE FIX
=======
Stop asking the classifier to be clever and make its input exact.  If the two
edges are BLOCK BOUNDARIES then a node on the leading edge has t = 0 exactly,
is never displaced, and the `any` rule lands on the true outline by itself --
there is nothing left to quantise.

Making them block boundaries does not need a new topology, a hybrid zone or an
embedded curve.  (Nor could it use one: gmsh honours `Curve In Surface` only
under the unstructured Delaunay and HXT algorithms, so an embedded curve and a
transfinite block are mutually exclusive.)  It needs the axial stations to stop
being planes:

    x  ->  x + d(x, r)

a shear in the (x, r) plane, applied to the finished node array exactly the way
finPatch applies the fin deformation.  Anchors, and the shift each one carries:

    fin_x_start()          0                     pinned
    FIN_LE_X_WALL          x_LE(r) - x_LE(R)     lands on the leading edge
    FIN_TE_X_WALL          x_TE(r) - x_TE(R)     lands on the trailing edge
    X_BASE                 relief * the above    keeps `tail` from collapsing
    X_BASE + X_WAKE_1      0                     pinned

d is linear in x between consecutive anchors, so inside the chord block

    x_new = (1 - phi) x_LE(r) + phi x_TE(r)

i.e. the axial index IS the chord fraction, at every span station.  60 cells
over 300 mm at the root and over 150 mm at the tip, with no extra machinery.

WHY THE WALL DOES NOT MOVE
==========================
Every shift is measured from its own value at r = R_BODY and G.fin_edge_x()
freezes below R_BODY, so d vanishes identically on the body surface -- on the
cylinder, on the boattail, and on the nose, which is outside the anchor range
anyway.  snap() has already put those nodes on the analytic surface to the
nanometre and this step leaves them there.

WHAT IT COSTS
=============
Skew, and honestly: a face lying on the leading edge is at 57.4 deg to the
axial direction, because that is the sweep.  No choice of blend improves it
while the radial lines are circles -- the only way down is a collar wrapped
around the planform edge, which is a different mesh generator.  Everything
else stays mild: the radial blend out to FIN_EDGE_R_BLEND contributes about
8 deg, and the four anchor intervals stretch or squeeze by

    cylfin    +411 %      finchord   -50 %      tail   -40 %      wake1  -17 %

all of which check() verifies are above the -100 % that would invert a cell.
"""

import numpy as np

import aconcaguaGeom as G
import meshParams as MP


def window(r):
    """Radial taper of the shear: 1 out to the fin tip, 0 by FIN_EDGE_R_BLEND.

    Smoothstep rather than linear so the shear has no kink at either end; the
    price is a peak slope 1.5x the average, which is folded into the ~8 deg
    quoted above.
    """
    r0 = G.FIN_TIP_R + max(MP.FIN_TIP_SMEAR, 0.0)
    r1 = MP.fin_edge_r_blend()
    t = np.clip((np.asarray(r, float) - r0) / (r1 - r0), 0.0, 1.0)
    return 1.0 - t * t * (3.0 - 2.0 * t)


def anchors():
    """The stations the shear moves, as (x_undeformed, edge, factor).

    `edge` is 'le', 'te' or None (pinned).  Ordered, and the first and last are
    always pinned so the shear is confined to the fin's own axial band.
    """
    return [(MP.fin_x_start(),           None, 0.0),
            (G.FIN_LE_X_WALL,            'le', 1.0),
            (G.FIN_TE_X_WALL,            'te', 1.0),
            (G.X_BASE,                   'te', float(MP.FIN_EDGE_TAIL_RELIEF)),
            (G.X_BASE + MP.X_WAKE_1,     None, 0.0)]


def shifts(r):
    """Shift of every anchor station at radius r.  Shape (n_anchor, len(r))."""
    le, te = G.fin_edge_x(r)
    w = window(r)
    s = {'le': (le - G.FIN_LE_X_WALL) * w,
         'te': (te - G.FIN_TE_X_WALL) * w,
         None: np.zeros_like(w)}
    return np.array([f * s[edge] for _, edge, f in anchors()])


def warp(nodes):
    """Shear x so the fin edges land on block boundaries.

    Returns (new_nodes, n_moved, max_shift).  Reads UNDEFORMED coordinates and
    must run BEFORE finPatch.deform(), which reads (x, r) and would otherwise
    be evaluating the planform at the wrong station.
    """
    xa = [a[0] for a in anchors()]
    p = nodes.copy()
    x = p[:, 0]
    r = np.hypot(p[:, 1], p[:, 2])
    band = (x >= xa[0]) & (x <= xa[-1])
    if not np.any(band):
        return p, 0, 0.0
    sa = shifts(r)
    d = np.zeros(len(p))
    for k in range(len(xa) - 1):
        m = band & (x >= xa[k]) & (x <= xa[k + 1])
        if not np.any(m):
            continue
        t = (x[m] - xa[k]) / (xa[k + 1] - xa[k])
        d[m] = sa[k][m] * (1.0 - t) + sa[k + 1][m] * t
    p[:, 0] = x + d
    return p, int((np.abs(d) > 1e-12).sum()), float(np.abs(d).max())


def check():
    """Assert the shear cannot invert a cell, and report how hard it pulls.

    An anchor interval of length L whose two shifts differ by ds is stretched
    by the factor 1 + ds / L.  Below zero it has folded the block over itself,
    and gmsh will not have warned because the mesh was valid when it was
    written -- the fold happens here, afterwards.  So check it here.
    """
    r = np.linspace(0.0, MP.fin_edge_r_blend() * 1.05, 4001)
    xa = [a[0] for a in anchors()]
    sa = shifts(r)
    worst = np.inf
    for k in range(len(xa) - 1):
        L = xa[k + 1] - xa[k]
        worst = min(worst, 1.0 + float(((sa[k + 1] - sa[k]) / L).min()))
    if worst <= 0.0:
        raise AssertionError(
            'FIN_EDGE_FIT would fold a block over itself (stretch factor '
            f'{worst:.3f} <= 0).  Raise FIN_EDGE_TAIL_RELIEF or FIN_X_LEAD.')
    return dict(min_stretch=worst)
