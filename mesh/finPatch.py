#!/usr/bin/env python3
"""
finPatch.py -- introduce the fins into the structured quadrant mesh.

Two steps, neither of which touches the block topology.

1. DEFORM the azimuthal coordinate of every node:

       theta -> theta_f(x,r) + theta * (90 - 2 theta_f) / 90

   A node that was on the symmetry plane (theta = 0) lands at theta_f, i.e. at
   z = -r sin(theta_f) = -t_half(x,r), which IS the fin surface.  The geometry
   comes out exact rather than approximated, the mesh stays conforming, and not
   one block is added.  It works only because the real fin is BEVELLED: t goes
   to zero continuously at the leading and trailing edges, so the deformation
   relaxes to nothing there.  The tip is the one genuine discontinuity and is
   smeared over FIN_TIP_SMEAR.

   Both symmetry planes get it, because the quadrant contains two half fins --
   one lying in each plane.  The full-thickness fin appears when two quadrants
   are stitched together (sectorAssembly.py): the two half-fins on either side
   of the interface plane are its two faces.

2. CLASSIFY the symmetry-plane faces.  A face whose nodes were all left in the
   plane is still symmetry; a face with any displaced node lies on the fin and
   becomes a wall.  The same rule decides, at assembly, which interface faces
   are stitched (merged into the interior) and which stay as fin walls -- the
   two decisions cannot disagree because they are the same test.
"""

import numpy as np

import aconcaguaGeom as G
import meshParams as MP


def fin_node_mask(nodes):
    """True for nodes the deformation will move: on the fin planform.

    Evaluated on UNDEFORMED coordinates.  Independent of which plane the node
    is on, so pass a plane mask and AND it in if you need one plane only.
    """
    section, smear = MP.FIN_SECTION, MP.FIN_TIP_SMEAR
    x, y, z = nodes[:, 0], nodes[:, 1], nodes[:, 2]
    r = np.hypot(y, z)
    out = np.zeros(len(nodes), bool)
    live = (r > 1e-9) & (x >= G.FIN_ROOT_LE - 1e-6) & (x <= G.FIN_TIP_TE + 1e-6)
    if np.any(live):
        t = G.fin_half_thickness(x[live], r[live], section, smear)
        out[live] = np.atleast_1d(t) > 1e-9
    return out


def deform(nodes):
    """Apply the fin deformation.  Returns (new_nodes, n_moved, max_displacement)."""
    section, smear = MP.FIN_SECTION, MP.FIN_TIP_SMEAR
    p = nodes.copy()
    x, y, z = p[:, 0], p[:, 1], p[:, 2]
    r = np.hypot(y, z)
    live = (r > 1e-9) & (x >= G.FIN_ROOT_LE - 1e-6) & (x <= G.FIN_TIP_TE + 1e-6)
    if not np.any(live):
        return p, 0, 0.0
    th = np.arctan2(-z[live], y[live])
    tf = np.atleast_1d(G.fin_half_angle(x[live], r[live], section, smear))
    half = np.pi / 2.0
    thn = tf + th * (half - 2.0 * tf) / half
    rr = r[live]
    ny, nz = rr * np.cos(thn), -rr * np.sin(thn)
    d = np.hypot(ny - y[live], nz - z[live])
    y[live], z[live] = ny, nz
    return p, int((d > 1e-12).sum()), float(d.max())


def split_symm(quads_symm, fin_nodes):
    """Split the symmetry-plane quads into (still symmetry, now fin)."""
    on_fin = fin_nodes[quads_symm].any(axis=1)
    return quads_symm[~on_fin], quads_symm[on_fin]


def wetted_area(nodes, quads):
    """Area of a quad patch, for checking the fin against its analytic value."""
    p = nodes[quads]
    d1 = p[:, 2] - p[:, 0]
    d2 = p[:, 3] - p[:, 1]
    return float(0.5 * np.linalg.norm(np.cross(d1, d2), axis=1).sum())


def wetted_area_analytic(n=2001):
    """True wetted area of ONE side of one fin: the area of the surface the
    deformation actually creates, z = -t_half(x, r).

    The reference the mesh converges to, and the one to compare against.  It
    exceeds the nominal planform for two reasons that are geometry, not
    discretisation:

      * the tip taper extends the fin FIN_TIP_SMEAR past FIN_TIP_R
      * the bevels and the tip taper are inclined, so the surface is larger
        than its projection by sqrt(1 + (dt/dx)^2 + (dt/dr)^2)

    For the as-drawn wedge fin that is 376.0 cm2 against a 361.5 cm2 nominal
    planform, a 4 % difference -- the same order as the discretisation error
    it was being blamed for.
    """
    section, smear = MP.FIN_SECTION, MP.FIN_TIP_SMEAR
    hi = G.FIN_TIP_R + max(smear, 0.0)
    x = np.linspace(G.FIN_ROOT_LE, G.FIN_TIP_TE, n)
    r = np.linspace(G.FIN_ROOT_R, hi, n)
    X, R = np.meshgrid(x, r, indexing='ij')

    def t(xx, rr):
        return np.atleast_1d(G.fin_half_thickness(xx, rr, section, smear)).reshape(xx.shape)

    h = 1e-6
    tt = t(X, R)
    tx = (t(X + h, R) - t(X - h, R)) / (2.0 * h)
    tr = (t(X, R + h) - t(X, R - h)) / (2.0 * h)
    dA = (x[1] - x[0]) * (r[1] - r[0])
    return float((np.sqrt(1.0 + tx ** 2 + tr ** 2) * (tt > 1e-12)).sum() * dA)
