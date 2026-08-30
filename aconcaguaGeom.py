#!/usr/bin/env python3
"""
aconcaguaGeom.py -- the GEOMETRY OF RECORD for the Aconcagua rocket.

This file replaces constant/triSurface/Aconcagua.stl as the authoritative
description of the shape.  The STL is a 2096-facet sample of these formulas;
these formulas are the thing itself.

Provenance
----------
Every constant below was measured back out of Aconcagua.stl and then verified
against the analytic form.  The nosecone was identified as a Von Karman
(LV-Haack, C = 0) ogive by fitting, not assumed:

    x = 0.11030  ->  STL r = 0.021781   LV-Haack r = 0.021768
    x = 0.40000  ->  STL r = 0.053387   LV-Haack r = 0.053396
    x = 0.62223  ->  STL r = 0.068972   LV-Haack r = 0.068971

Agreement to ~5 significant figures over the whole nose is not a coincidence;
it is the curve the CAD was drawn with.  run `python3 aconcaguaGeom.py` to
re-run that verification against the STL at any time.

Why this matters for meshing
----------------------------
A triangulated surface has no tangent plane and no curvature.  Every mesher
that consumes one must *guess* both, which is why snappyHexMesh needs
surfaceFeatureExtract, resolveFeatureAngle, snap tolerances and nSolveIter --
those knobs exist to paper over information the STL threw away.  With r(x) in
closed form we differentiate instead of guessing: wall normals are exact,
curvature-based sizing is exact, and a node can be placed *on* the surface
rather than iterated towards it.
"""

import numpy as np

# ----------------------------------------------------------------- constants --
# THE FIVE NUMBERS THAT DEFINE THE BODY.  Everything else is DERIVED from them,
# so changing one cannot leave another stale.
L_NOSE   = 0.80000      # m  nose length
L_CYL    = 2.03000      # m  cylindrical section length
L_TAIL   = 0.12500      # m  boattail length
R_BODY   = 0.07550      # m  cylinder radius   (D = 0.15100)
R_BASE   = 0.05500      # m  base radius

# --- derived.  Do NOT edit: they follow. --------------------------------------
X_BODY_1 = L_NOSE                       # nose / cylinder junction
X_BODY_2 = X_BODY_1 + L_CYL             # cylinder / boattail junction
X_BASE   = X_BODY_2 + L_TAIL            # base plane == total length
L_TOTAL  = X_BASE
D_BODY   = 2.0 * R_BODY

# An earlier version carried L_NOSE and X_BODY_1 as two INDEPENDENT constants.
# Changing L_NOSE alone then opened an 8.90 mm step in r_body at the junction --
# 29 first-cell heights -- with no error and no warning, because each formula
# was individually correct.  Deriving the junction removes the failure mode
# instead of documenting it.  validate() below asserts what is left.

# --- fins: 4 planar clipped deltas at theta = 0, 90, 180, 270 -----------------
# Anchored to the BASE rather than to absolute x, so they travel with the tail
# when the body length changes.
N_FINS       = 8
FIN_T        = 0.012    # m  full thickness
FIN_ROOT_R   = 0.07500  # m  radial station of the root chord
FIN_TIP_R    = 0.23550  # m  semi-span
FIN_ROOT_LE  = X_BASE - 0.42578
FIN_ROOT_TE  = X_BASE - 0.12531
FIN_TIP_LE   = X_BASE - 0.17500
FIN_TIP_TE   = X_BASE - 0.02500
FIN_LE_BEVEL = 0.03320  # m  sharp LE -> full thickness (ABSOLUTE, not % chord)
FIN_TE_BEVEL = 0.02116  # m  full thickness -> sharp TE

# ------------------------------------------------------------ nose profile ---
def r_nose(x):
    """Von Karman / LV-Haack (C = 0) ogive radius, x in [0, L_NOSE].

        theta = arccos(1 - 2x/L)
        r     = R/sqrt(pi) * sqrt(theta - sin(2 theta)/2)

    Near the tip this behaves as r ~ x**(3/4): the slope is INFINITE at x = 0.
    That is a real property of the shape, not a numerical artefact, and it is
    the single fact that dictates the mesh topology -- see the note on the
    butterfly cap in buildHexBody.py.  A radial O-grid cannot be carried into
    a point of vertical tangency without collapsing.
    """
    x = np.asarray(x, dtype=float)
    t = np.arccos(np.clip(1.0 - 2.0 * x / L_NOSE, -1.0, 1.0))
    return R_BODY / np.sqrt(np.pi) * np.sqrt(np.maximum(t - 0.5 * np.sin(2.0 * t), 0.0))


def drdx_nose(x, h=1e-9):
    """d r / d x on the nose, by central difference on the closed form."""
    x = np.asarray(x, dtype=float)
    xs = np.clip(x, h, L_NOSE - h)
    return (r_nose(xs + h) - r_nose(xs - h)) / (2.0 * h)


def x_nose_of_r(r):
    """Inverse of r_nose: the x at which the nose has radius r.

    Monotone, so a bisection is exact to machine precision and needs no
    Newton safeguard against the infinite slope at the tip.
    """
    r = np.atleast_1d(np.asarray(r, dtype=float))
    lo = np.zeros_like(r)
    hi = np.full_like(r, L_NOSE)
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        go = r_nose(mid) < r
        lo = np.where(go, mid, lo)
        hi = np.where(go, hi, mid)
    out = 0.5 * (lo + hi)
    return out if out.size > 1 else float(out[0])


# ------------------------------------------------------------ full profile ---
def r_body(x):
    """Radius of the body of revolution at station x, for x in [0, X_BASE].

    Piecewise: LV-Haack ogive | cylinder | straight conical frustum.
    """
    x = np.atleast_1d(np.asarray(x, dtype=float))
    r = np.empty_like(x)

    nose = x <= X_BODY_1
    cyl  = (x > X_BODY_1) & (x <= X_BODY_2)
    tail = x > X_BODY_2

    r[nose] = r_nose(x[nose])
    r[cyl]  = R_BODY
    r[tail] = R_BODY + (R_BASE - R_BODY) * (x[tail] - X_BODY_2) / (X_BASE - X_BODY_2)
    return r if r.size > 1 else float(r[0])


def drdx_body(x):
    x = np.atleast_1d(np.asarray(x, dtype=float))
    s = np.empty_like(x)
    nose = x <= X_BODY_1
    cyl  = (x > X_BODY_1) & (x <= X_BODY_2)
    tail = x > X_BODY_2
    s[nose] = drdx_nose(x[nose])
    s[cyl]  = 0.0
    s[tail] = (R_BASE - R_BODY) / (X_BASE - X_BODY_2)
    return s if s.size > 1 else float(s[0])


def wall_normal(x):
    """Outward unit normal of the surface of revolution in the (x, r) plane."""
    s = np.asarray(drdx_body(x), dtype=float)
    n = np.sqrt(1.0 + s * s)
    return (-s / n, 1.0 / n)          # (n_x, n_r)


def x_at_slope(target):
    """Station on the NOSE where dr/dx equals `target` (slope decreases with x)."""
    lo, hi = 1e-9, L_NOSE
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if drdx_nose(mid) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


BOATTAIL_HALF_ANGLE = float(np.degrees(np.arctan2(R_BODY - R_BASE, X_BASE - X_BODY_2)))


# ------------------------------------------------------------- fin sections --
def fin_planform(sigma):
    """Leading edge x and chord at spanwise fraction sigma in [0, 1].

    Straight-tapered clipped delta, so both are linear in sigma.  sigma = 0 is
    the root chord at r = FIN_ROOT_R, sigma = 1 the tip at r = FIN_TIP_R.
    """
    sigma = np.asarray(sigma, dtype=float)
    le = FIN_ROOT_LE + sigma * (FIN_TIP_LE - FIN_ROOT_LE)
    te = FIN_ROOT_TE + sigma * (FIN_TIP_TE - FIN_ROOT_TE)
    return le, te - le


def fin_halfthickness(xi, chord, section='wedge', t=None, **kw):
    """Half-thickness of the fin section at chord fraction xi in [0, 1].

    `section` selects the profile.  The fin is the one place on this vehicle
    where the section is a genuine design choice rather than a measurement, so
    it is a parameter and not a constant:

      'wedge'    the as-drawn shape: sharp LE, straight ramp to full thickness
                 over FIN_LE_BEVEL, flat, straight ramp to a sharp TE over
                 FIN_TE_BEVEL.  Bevels are ABSOLUTE lengths, so the section
                 changes shape from root to tip -- 11 % of chord at the root,
                 22 % at the tip.  That is what a machined fin does.
      'diamond'  pure double wedge, maximum thickness at `xmax` (default 0.5)
      'biconvex' circular-arc, the classic supersonic section
      'naca'     NACA symmetric 4-digit, thickness ratio t/chord
      'naca_te'  same but with the TE closed exactly (last coefficient adjusted)

    Every one of these is analytic, which is the entire point: the surface, its
    tangent and its curvature are all available in closed form, so leading-edge
    cell size can be set FROM the leading-edge radius instead of guessed.
    """
    xi = np.clip(np.asarray(xi, dtype=float), 0.0, 1.0)
    tt = FIN_T if t is None else t
    h = 0.5 * tt

    if section == 'wedge':
        a = FIN_LE_BEVEL / chord
        b = 1.0 - FIN_TE_BEVEL / chord
        return h * np.clip(np.minimum(xi / a, (1.0 - xi) / (1.0 - b)), 0.0, 1.0)

    if section == 'diamond':
        m = kw.get('xmax', 0.5)
        return h * np.where(xi < m, xi / m, (1.0 - xi) / (1.0 - m))

    if section == 'biconvex':
        return h * 4.0 * xi * (1.0 - xi)

    if section in ('naca', 'naca_te'):
        c4 = -0.1036 if section == 'naca_te' else -0.1015   # closed vs open TE
        return (tt / 0.20) * chord * (0.29690 * np.sqrt(xi) - 0.12600 * xi
                                      - 0.35160 * xi ** 2 + 0.28430 * xi ** 3
                                      + c4 * xi ** 4)
    raise ValueError(f'unknown fin section {section!r}')


def fin_le_radius(chord, section='naca', t=None):
    """Leading-edge radius, which is what actually sets the LE cell size.

    A sharp LE ('wedge', 'diamond', 'biconvex') has radius 0 and needs an
    H- or O-grid closing on the edge.  A NACA LE is round and wants a C-grid
    with about 20 cells over the nose; r_LE = 1.1019 (t/c)^2 c.
    """
    if section not in ('naca', 'naca_te'):
        return 0.0
    tt = FIN_T if t is None else t
    return 1.1019 * (tt / chord) ** 2 * chord


def fin_half_thickness(x, r, section='wedge', tip_smear=0.004):
    """Fin half-thickness at body station x and radius r.  Zero off the planform.

    This is the field the mesh deformation reads.  It goes to zero CONTINUOUSLY
    at the leading and trailing edges, because the real fin is bevelled there --
    which is what lets the fin be introduced without a topology change.  Only
    the tip is a genuine discontinuity, so it is smeared over `tip_smear`.
    """
    x = np.atleast_1d(np.asarray(x, float))
    r = np.atleast_1d(np.asarray(r, float))
    out = np.zeros(np.broadcast(x, r).shape)

    span = FIN_TIP_R - FIN_ROOT_R
    sigma = np.clip((r - FIN_ROOT_R) / span, 0.0, 1.0)
    le, chord = fin_planform(sigma)
    xi = (x - le) / chord
    xb = np.broadcast_to(x, out.shape)
    rb = np.broadcast_to(r, out.shape)
    xib = np.broadcast_to(xi, out.shape)
    cb = np.broadcast_to(chord, out.shape)
    hi = FIN_TIP_R + max(tip_smear, 0.0)
    # no lower radial bound: below the root chord the fin is buried in the body,
    # and the STL models it that way (full-thickness vertices sit at r = 0.07494)
    inside = (xib >= 0.0) & (xib <= 1.0) & (rb <= hi + 1e-12)
    if np.any(inside):
        out[inside] = np.array([fin_halfthickness(float(a), float(b), section)
                                for a, b in zip(xib[inside], cb[inside])])
        if tip_smear > 0.0:      # taper the tip instead of ending on a cliff
            out *= np.clip((hi - rb) / tip_smear, 0.0, 1.0)
    return out if out.size > 1 else float(out.reshape(-1)[0])


def fin_half_angle(x, r, section='wedge', tip_smear=0.004):
    """Azimuthal half-angle the fin subtends at (x, r).  arcsin(t_half / r)."""
    t = fin_half_thickness(x, r, section, tip_smear)
    r = np.maximum(np.asarray(r, float), 1e-9)
    return np.arcsin(np.clip(t / r, 0.0, 0.999))


def on_fin_planform(x, r, section='wedge'):
    """True where a point of the symmetry plane is covered by the fin."""
    return fin_half_thickness(x, r, section, tip_smear=0.0) > 1e-9


# --------------------------------------------------------------- self-check --
def validate():
    """Assert the invariants the mesh topology depends on.  Cheap insurance."""
    errs = []
    if abs(r_nose(L_NOSE) - R_BODY) > 1e-12:
        errs.append('nose does not close on the cylinder radius')
    xs = np.linspace(1e-9, L_NOSE, 20000)
    if np.any(np.diff(r_nose(xs)) < 0):
        errs.append('r_nose is not monotone -- the butterfly cap map is not a bijection')
    if not (0.0 < FIN_ROOT_LE < FIN_ROOT_TE <= X_BASE):
        errs.append('fin root chord falls outside the body')
    if FIN_ROOT_R > r_body(FIN_ROOT_TE) + 1e-9:
        errs.append('fin root radius is outside the body surface it mounts on')
    if FIN_LE_BEVEL + FIN_TE_BEVEL >= (FIN_TIP_TE - FIN_TIP_LE):
        errs.append('fin bevels overlap at the tip chord')
    if errs:
        raise AssertionError('geometry is inconsistent:\n  - ' + '\n  - '.join(errs))
    return True


# --------------------------------------------------------- reference values --
def reference_values(quarter=True):
    """Aref / lRef for force coefficients.  Quarter symmetry divides Aref by 4."""
    A = np.pi * D_BODY ** 2 / 4.0
    return dict(Aref=A / (4.0 if quarter else 1.0), lRef=D_BODY, Sref_full=A)


# ------------------------------------------------------------ verification ---
def verify_against_stl(path='constant/triSurface/Aconcagua.stl', verbose=True):
    """Re-derive the fit from the STL.  Returns max |r_stl - r_analytic|."""
    import re
    txt = open(path).read()
    solids = dict((n, np.array([[float(c) for c in m] for m in
                   re.findall(r'vertex\s+(\S+)\s+(\S+)\s+(\S+)', b)]))
                  for n, b in re.findall(r'solid (\S+)(.*?)endsolid', txt, re.S))

    worst = 0.0
    rows = []
    for name in ('nosecone', 'body', 'boattail'):
        v = solids[name]
        x = v[:, 0]
        r = np.hypot(v[:, 1], v[:, 2])
        ra = r_body(x)
        e = np.abs(r - ra)
        worst = max(worst, e.max())
        rows.append((name, len(v) // 3, e.max(), e.max() / R_BODY))

    if verbose:
        print(f'{"solid":<10}{"facets":>8}{"max |dr| [m]":>16}{"rel. to R":>14}')
        for n, f, e, rel in rows:
            print(f'{n:<10}{f:>8}{e:>16.3e}{rel:>14.2e}')
        print(f'\nworst deviation over the whole body : {worst:.3e} m '
              f'({worst / R_BODY * 100:.4f} % of body radius)')
        print(f'boattail half-angle                 : {BOATTAIL_HALF_ANGLE:.3f} deg')
        print(f'slope dr/dx = 1 (45 deg) at x       : {x_at_slope(1.0):.6f} m, '
              f'r = {r_body(x_at_slope(1.0)):.6f} m')
    return worst


if __name__ == '__main__':
    import sys
    validate()
    p = sys.argv[1] if len(sys.argv) > 1 else 'constant/triSurface/Aconcagua.stl'
    verify_against_stl(p)
