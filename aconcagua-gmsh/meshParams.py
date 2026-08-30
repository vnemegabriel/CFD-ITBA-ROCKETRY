#!/usr/bin/env python3
"""
meshParams.py -- the ONLY file you normally edit.

The interface here is deliberately different from a snappyHexMeshDict.  You do
not specify refinement LEVELS and discover the cell size afterwards; you
specify the CELL SIZE at each end of each mesh segment, and the cell count
falls out of the geometry:

    h_N = h_1 c^(n-1)                 (geometric stack)
    L   = h_1 (c^n - 1)/(c - 1)
    =>  c = (L - h_1)/(L - h_N),   n = 1 + ln(h_N/h_1)/ln(c)

Both equations are closed form, so "0.3 mm at the wall growing to 20 mm at
r = 0.35 m" is a statement about the mesh, not a wish about it.  Nothing here
is iterated, relaxed, or subject to a quality-driven retreat.
"""

import math
import os

import numpy as np
import aconcaguaGeom as G

# ====================== GLOBAL COARSITY -- ONE NUMBER ========================
# Multiplies every cell SIZE in this file and divides every cell COUNT, so the
# whole mesh scales consistently and every distribution re-solves from the new
# sizes.  Cells go roughly as 1/H_SCALE^3.
#
#   H_SCALE = 1.0   as tabulated below
#   H_SCALE = 2.0   about 1/8 the cells, for a quick look
#   H_SCALE = 0.7   about 3x the cells, for a convergence study
#
# H_SCALE_WALL = False holds y1 (and therefore y+) fixed while everything else
# scales -- what you want for a grid-convergence study with wall functions.
H_SCALE      = 1.0
H_SCALE_WALL = True

H_SCALE = float(os.environ.get('AG_COARSE', H_SCALE))    # env override for smoke tests


def _h(v):
    """Scale a cell size."""
    return None if v is None else v * H_SCALE


def _n(v, lo=1):
    """Scale a cell count."""
    return max(lo, int(round(v / H_SCALE)))

# ------------------------------------------------------------------- flow ----
U, NU, RHO = 100.0, 1.5e-5, 1.225
YPLUS_TARGET = 32.0          # first cell CENTRE; wall-function range

# -------------------------------------------------------------- azimuthal ---
# The quadrant is divided into N_AZ_BLOCKS azimuthal BLOCKS, each carrying
# N_AZ_CELLS cells.  Block WIDTHS are graded toward both symmetry planes --
# that is where the fins are -- while each block stays internally uniform.
#
# Grading by block width rather than by a progression across the whole 90 deg
# is the point.  A transfinite quad blends the radius LINEARLY IN INDEX, so a
# clustered distribution decouples angle from index and the interpolation error
# grows with the angular span of the block.  Keep the spans small and the error
# is bounded by the block's own arc-to-chord sagitta, which is tiny.
N_AZ_BLOCKS     = 12         # per quadrant; MUST be even (the butterfly core
                             # is an (N/2) x (N/2) grid of sub-blocks)
N_AZ_CELLS      = 3          # cells per block -> 36 per quadrant, 144 around
AZ_BLOCK_GROWTH = 1.50       # width ratio between successive blocks, from the
                             # symmetry planes inward
AZ_FIN_H        = 0.5e-3     # m, first azimuthal cell AT r = R_BODY inside the
                             # block touching each symmetry plane.  None = that
                             # block is uniform too.
N_RING  = 10                 # cells across the butterfly ring


def n_az_cells():
    return _n(N_AZ_CELLS, 1)


def n_ring():
    return _n(N_RING, 2)

# ------------------------------------------------------------------- fins ---
# The fin is introduced by DEFORMING the azimuthal coordinate, not by adding
# blocks: theta -> theta_f + theta (90 - 2 theta_f)/90, with theta_f(x,r) =
# arcsin(t_half/r).  Nodes at theta = 0 land exactly on z = -t_half, which IS
# the fin surface, so the geometry is exact and the topology is untouched.
# It works only because the real fin is bevelled: t -> 0 continuously at the
# leading and trailing edges, so the deformation relaxes to zero there.
FINS_ON      = True
FIN_SECTION  = 'wedge'       # wedge | diamond | biconvex | naca | naca_te
FIN_TIP_SMEAR = 0.004        # m, radial band over which the tip taper closes

# --- refinement around the fins ---------------------------------------------
# Three directions, three knobs.
#   azimuthal  AZ_FIN_H        first cell off the fin surface
#              AZ_BLOCK_GROWTH how fast blocks widen away from it
#   streamwise FIN_H_X         cell size over the fin chord; splits `cyl` so the
#                              refined block starts FIN_X_LEAD ahead of the root
#                              leading edge and runs to the base.  None = off.
#   radial     add a ZONE_R entry just outside FIN_TIP_R (0.2355 m) with a small
#              ZONE_H -- see the zone block above.  Adding a zone also means
#              setting WAKE_ZONE_K to whichever zone should open into the wake.
FIN_H_X      = 0.005         # m   None = leave `cyl` and `tail` unrefined
FIN_X_LEAD   = 0.060         # m   of cylinder ahead of the root LE to include

# ======================= REFINEMENT ZONES -- EDIT HERE =======================
# Concentric zones outward from the wall.  ZONE_R[k] is where zone k ENDS, so
# the list must increase and the LAST ENTRY IS THE FARFIELD RADIUS.  There is
# no second place to state that: R_FAR below is derived from this list, not a
# parameter.  (It used to be an independent constant that nothing read, so
# editing it changed the documentation and not the mesh.)
#
# ZONE_H[k] is the cell size at the OUTER edge of zone k.  The inner edge of
# zone 0 is y1, from the y+ target.
ZONE_R = [0.35, 2.00, 11.775]      # m
ZONE_H = [0.020, 0.200, 1.600]     # m

# AMPLIFIED WAKE.  Zone 0 is the fine radial band; downstream of the base it
# OPENS OUT so the band tracks the spreading wake instead of the wake growing
# out of it.  A turbulent wake spreads as x^(1/2), hence the exponent.
WAKE_ZONE_K   = 0            # index of the zone that opens into the wake
ZONE0_R_WAKE  = 0.70         # m   that zone's outer radius at the outlet
WAKE_SPREAD_P = 0.5

# Internal blend tubes at the two open ends, as FRACTIONS OF ZONE 0.  They are
# fractions and not metres so they cannot outgrow the zone that contains them --
# a failure mode that used to produce inverted blocks with no warning.
F_INLET     = 0.63           # of ZONE_R[0]          -> 0.221 m
F_WAKE_OUT  = 0.50           # of the zone-0 radius at the outlet -> 0.350 m

# --- derived: do NOT edit ----------------------------------------------------
R_FAR  = ZONE_R[-1]
SHELLS = [dict(r_out=r, h_in=(None if k == 0 else ZONE_H[k - 1]), h_out=ZONE_H[k])
          for k, r in enumerate(ZONE_R)]

# ------------------------------------------------------------- streamwise ---
# h_start / h_end in metres; None means "continue from the neighbouring
# segment", so you cannot accidentally create a jump at a block interface.
SEGMENTS = [
    ('up',    dict(h_start=None,   h_end=1.5e-3)),    # inlet -> butterfly cap
    ('nose',  dict(h_start=1.5e-3, h_end=0.012)),     # cap   -> x = 0.80
    ('cyl',   dict(h_start=0.012,  h_end=0.012)),     # cylinder
    ('tail',  dict(h_start=0.008,  h_end=0.008)),     # boattail
    ('wake1', dict(h_start=None,   h_end=0.008)),     # base  -> +0.30 m
    ('wake2', dict(h_start=0.008,  h_end=0.030)),     # +0.30 -> +4.00 m
    ('wake3', dict(h_start=0.030,  h_end=1.500)),     # +4.00 -> outlet
]
# The inlet plane has no neighbouring segment to inherit a size from.  Given as
# a FRACTION of the first upstream sub-block rather than in metres, so a global
# coarsening cannot ask for a cell larger than the block it has to fit inside.
F_UP_INLET   = 0.50          # of (x_cap - x_in) / N_UP_BLOCKS  -> 1.48 m at H_SCALE 1
F_UP_MAX     = 0.60          # hard ceiling on that fraction
H_WAKE_BASE  = 5.0e-4        # m  first streamwise cell off the flat base

# --------------------------------------------------------- butterfly cap ----
# The LV-Haack tip has INFINITE dr/dx, so a radial O-grid cannot be carried
# into it -- the topology must change.  It changes here, and nowhere else.
CAP_R_FRAC   = 0.10          # cap rim radius / R_BODY
CORE_FRAC    = 0.45          # core square half-width / local inner radius
# The upstream inner boundary leaves the cap rim TANGENT to the wall:
#     Ri(x)^2 = r_cap^2 + 2 r_cap s_cap (x_cap - x) + lambda (x_cap - x)^2
# with lambda fixed by Ri(x_in) = R_INLET.  A simple power blend r_cap +
# dR * xi^q was tried first and gave 77 deg non-orthogonality at the rim,
# because xi^q has UNBOUNDED CURVATURE at xi = 0 for 1 < q < 2 -- the surface
# is tangent-discontinuous in the second derivative exactly where it matters.

# ------------------------------------------------------------ domain size ---
UPSTREAM_L, DOWNSTREAM_L = 6.0, 13.0        # body lengths; matches the snappy box
X_WAKE_1     = 0.30                          # m behind the base
X_WAKE_2     = 4.00          # near wake now runs to 27 body diameters

QUARTER = True
PATCHES = dict(inlet='inlet', outlet='outlet', farfield='box', symmetry='symm',
               nose='cone', body='walls', tail='tail')


# ------------------------------------------------------------------ solver --
def geometric(L, h1, hN):
    """Cell count and ratio for a geometric stack of length L from h1 to hN."""
    if h1 is None or hN is None:
        raise ValueError('both end sizes must be resolved before solving')
    if h1 >= L or hN >= L:
        raise ValueError(f'cell size ({h1:g}, {hN:g}) is not smaller than the '
                         f'segment it has to fill ({L:g} m)')
    if abs(h1 - hN) / max(h1, hN) < 1e-9:
        return max(int(round(L / h1)), 1), 1.0
    c = (L - h1) / (L - hN)
    n = max(int(round(1.0 + np.log(hN / h1) / np.log(c))), 2)
    return n, refit(L, n, h1)


def stack_sum(c, n):
    """(c^n - 1)/(c - 1), the length of a geometric stack of unit first cell.

    Overflow-safe.  A direct `c ** n` raises OverflowError once n log c passes
    ~709, which is not an exotic case: asking for y+ = 1 with a gentle
    expansion ratio drives n into the hundreds and the solver used to die with
    a bare `OverflowError: numerical result out of range` from inside a lambda,
    with nothing to say which parameter caused it.  Returning inf keeps the
    bisection monotone and lets it converge from the other side.
    """
    if abs(c - 1.0) < 1e-12:
        return float(n)
    ln = n * math.log(c)
    if ln > 700.0:
        return math.inf
    return (math.exp(ln) - 1.0) / (c - 1.0)


def refit(L, n, h1):
    """Ratio giving exactly n cells of total length L starting at h1."""
    target = L / h1
    if abs(stack_sum(1.0, n) - target) < 1e-9:
        return 1.0
    a, b = (1.0, 50.0) if target > n else (1e-4, 1.0)
    fa = stack_sum(a, n) - target
    for _ in range(300):
        m = 0.5 * (a + b)
        fm = stack_sum(m, n) - target
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    return 0.5 * (a + b)


def derived():
    validate_params()
    L = G.L_TOTAL
    Re_L = U * L / NU
    Cf   = 0.0576 * Re_L ** -0.2
    utau = U * np.sqrt(Cf / 2.0)
    y1   = 2.0 * YPLUS_TARGET * NU / utau        # first cell HEIGHT
    if H_SCALE_WALL:
        y1 *= H_SCALE
    delta = 0.37 * L * Re_L ** -0.2

    r_cap = CAP_R_FRAC * G.R_BODY
    x_cap = G.x_nose_of_r(r_cap)

    x_in  = -UPSTREAM_L * L
    x_out = G.X_BASE + DOWNSTREAM_L * L

    # radial shells, resolved at the CYLINDER station where the wall is
    shells, h_prev, r_prev = [], y1, G.R_BODY
    for s in SHELLS:
        h_in = _h(s['h_in']) if s['h_in'] is not None else h_prev
        h_out = _h(s['h_out'])
        n, c = geometric(s['r_out'] - r_prev, h_in, h_out)
        shells.append(dict(r_out=s['r_out'], n=n, c=c, h_in=h_in, h_out=h_out))
        h_prev, r_prev = h_out, s['r_out']

    # cells inside the boundary layer
    cum, n_delta, h = 0.0, 0, y1
    for _ in range(shells[0]['n']):
        cum += h
        if cum <= delta:
            n_delta += 1
        h *= shells[0]['c']

    return dict(Re_L=Re_L, Cf=Cf, utau=utau, y1=y1, delta=delta,
                x_in=x_in, x_out=x_out, r_cap=r_cap, x_cap=x_cap,
                cap_angle=float(np.degrees(np.arctan(G.drdx_body(x_cap)))),
                shells=shells, n_rad=sum(s['n'] for s in shells),
                n_delta=n_delta,
                ds=2.0 * np.pi * G.R_BODY / (4.0 * N_AZ_BLOCKS * n_az_cells()),
                n_circ=4 * N_AZ_BLOCKS * n_az_cells())


def fin_x_start():
    """Where the streamwise fin-refined block begins."""
    return G.FIN_ROOT_LE - FIN_X_LEAD


def segment_bounds(d):
    """Streamwise segment end stations, in order."""
    XB = G.X_BASE
    out = [('up', d['x_in'], d['x_cap']), ('nose', d['x_cap'], G.X_BODY_1)]
    if FIN_H_X is not None:
        out += [('cyl', G.X_BODY_1, fin_x_start()),
                ('cylfin', fin_x_start(), G.X_BODY_2)]
    else:
        out += [('cyl', G.X_BODY_1, G.X_BODY_2)]
    out += [('tail',  G.X_BODY_2, XB),
            ('wake1', XB,         XB + X_WAKE_1),
            ('wake2', XB + X_WAKE_1, XB + X_WAKE_2),
            ('wake3', XB + X_WAKE_2, d['x_out'])]
    return out


def h_up_inlet(d):
    """Streamwise cell at the inlet plane, in metres.

    A fraction of the first upstream sub-block, scaled by H_SCALE and capped, so
    it stays inside the block whatever the global coarsity is set to.
    """
    n = max(2, int(round(N_UP_BLOCKS / math.sqrt(H_SCALE))))
    L = (d['x_cap'] - d['x_in']) / n
    return min(F_UP_INLET * H_SCALE, F_UP_MAX) * L


def segment_sizes():
    """Cell sizes per streamwise segment, with the fin overrides folded in."""
    spec = dict(SEGMENTS)
    spec['wake1'] = dict(spec['wake1'], h_start=H_WAKE_BASE)
    if FIN_H_X is not None:
        spec['cyl'] = dict(spec['cyl'], h_end=FIN_H_X)
        spec['cylfin'] = dict(h_start=FIN_H_X, h_end=FIN_H_X)
        spec['tail'] = dict(h_start=FIN_H_X, h_end=FIN_H_X)
    return {k: dict(h_start=_h(v['h_start']), h_end=_h(v['h_end']))
            for k, v in spec.items()}


def axial_plan(d):
    """Resolve every streamwise segment to (n, ratio), honouring the Nones."""
    spec = segment_sizes()
    spec['up'] = dict(spec['up'], h_start=h_up_inlet(d))
    out = {}
    for name, x0, x1 in segment_bounds(d):
        n, c = geometric(x1 - x0, spec[name]['h_start'], spec[name]['h_end'])
        out[name] = dict(n=n, c=c, x0=x0, x1=x1,
                         h0=spec[name]['h_start'], h1=spec[name]['h_end'])
    return out


def zone_r_out(k, x):
    """Outer radius of zone k at station x.

    Zone WAKE_ZONE_K opens out downstream of the base so the fine radial band
    tracks the spreading wake; the others are cylinders.
    """
    if k != WAKE_ZONE_K or x <= G.X_BASE:
        return ZONE_R[k]
    x_out = G.X_BASE + DOWNSTREAM_L * G.L_TOTAL
    xi = min(max((x - G.X_BASE) / (x_out - G.X_BASE), 0.0), 1.0)
    return ZONE_R[k] + (ZONE0_R_WAKE - ZONE_R[k]) * xi ** WAKE_SPREAD_P


def r_inlet():
    """Blend-tube radius at the inlet plane.  A fraction of zone 0, so it
    cannot outgrow the zone it lives inside."""
    return F_INLET * ZONE_R[0]


def r_wake_out():
    """Blend-tube radius at the outlet plane, likewise a fraction."""
    return F_WAKE_OUT * zone_r_out(0, G.X_BASE + DOWNSTREAM_L * G.L_TOTAL)


# ---------------------------------------------------------------- azimuthal --
def az_angles():
    """Azimuthal BLOCK boundaries, 0 .. 90 deg, in radians.

    Widths grow from both symmetry planes inward by AZ_BLOCK_GROWTH, mirrored
    about 45 deg, so the finest blocks sit where the fins are.
    """
    n = N_AZ_BLOCKS
    if n % 2:
        raise ValueError(f'N_AZ_BLOCKS must be even, got {n}')
    half = n // 2
    q = float(AZ_BLOCK_GROWTH)
    w = np.array([q ** k for k in range(half)], dtype=float)
    w = np.concatenate([w, w[::-1]])
    w *= (np.pi / 2.0) / w.sum()
    return np.concatenate([[0.0], np.cumsum(w)])


def az_coefs():
    """Progression coefficient for each azimuthal block, in canonical
    (increasing theta) direction.  Only the two blocks touching a symmetry
    plane are graded; the rest are uniform."""
    th = az_angles()
    c = [1.0] * N_AZ_BLOCKS
    nc = n_az_cells()
    if AZ_FIN_H is not None:
        arc0 = G.R_BODY * (th[1] - th[0])
        if _h(AZ_FIN_H) < arc0 / nc:
            k = refit(arc0, nc, _h(AZ_FIN_H))
            c[0] = k                       # cluster at theta = 0
            c[-1] = 1.0 / k                # cluster at theta = 90
    return c


# ------------------------------------------------------------- consistency --
def validate_params():
    """Check the couplings between zone radii BEFORE anything is built.

    Every one of these used to be either silent or a confusing error from deep
    inside the cell-size solver.  They are cheap; run them always.
    """
    e = []
    if len(ZONE_R) != len(ZONE_H):
        e.append(f'ZONE_R has {len(ZONE_R)} entries, ZONE_H has {len(ZONE_H)}')
    if any(b <= a for a, b in zip(ZONE_R[:-1], ZONE_R[1:])):
        e.append(f'ZONE_R must increase: {ZONE_R}')
    if ZONE_R and ZONE_R[0] <= G.R_BODY:
        e.append(f'ZONE_R[0] = {ZONE_R[0]} is inside the body (R_BODY = {G.R_BODY})')
    k = WAKE_ZONE_K
    if not 0 <= k < len(ZONE_R):
        e.append(f'WAKE_ZONE_K = {k} is not a zone index')
    else:
        if ZONE0_R_WAKE < ZONE_R[k]:
            e.append(f'ZONE0_R_WAKE {ZONE0_R_WAKE} < ZONE_R[{k}] {ZONE_R[k]}: '
                     f'zone {k} would close up downstream instead of opening out')
        if k + 1 < len(ZONE_R) and ZONE0_R_WAKE >= ZONE_R[k + 1]:
            e.append(f'ZONE0_R_WAKE {ZONE0_R_WAKE} >= ZONE_R[{k+1}] {ZONE_R[k+1]}: '
                     f'zone {k} would swallow zone {k+1} at the outlet. Either lower '
                     f'ZONE0_R_WAKE or set WAKE_ZONE_K to the outermost fine zone.')
    if FIN_H_X is not None and not G.X_BODY_1 < fin_x_start() < G.FIN_ROOT_LE:
        e.append(f'FIN_X_LEAD = {FIN_X_LEAD} puts the fin block start at '
                 f'{fin_x_start():.4f}, outside the cylinder '
                 f'({G.X_BODY_1} .. {G.FIN_ROOT_LE})')
    if H_SCALE <= 0:
        e.append(f'H_SCALE must be > 0, got {H_SCALE}')
    if not 0.0 < F_INLET < 1.0:
        e.append(f'F_INLET must be in (0,1), got {F_INLET}')
    if not 0.0 < F_WAKE_OUT < 1.0:
        e.append(f'F_WAKE_OUT must be in (0,1), got {F_WAKE_OUT}')
    if not 0.0 < F_UP_INLET < 1.0 or not 0.0 < F_UP_MAX < 1.0:
        e.append('F_UP_INLET and F_UP_MAX must be in (0,1)')
    if not 0.0 < CORE_FRAC < 1.0 / np.sqrt(2.0):
        e.append(f'CORE_FRAC must be < 1/sqrt(2) = 0.707 or the core square '
                 f'pokes through its own ring; got {CORE_FRAC}')
    if N_AZ_BLOCKS % 2 or N_AZ_BLOCKS < 2:
        e.append(f'N_AZ_BLOCKS must be even and >= 2, got {N_AZ_BLOCKS}')
    for k, b in enumerate(ZONE_H):
        span = ZONE_R[k] - (G.R_BODY if k == 0 else ZONE_R[k - 1])
        if _h(b) >= span:
            e.append(f'ZONE_H[{k}] = {b} m (x H_SCALE = {_h(b):.4g}) is not smaller '
                     f'than zone {k}, which is {span:.4f} m wide')
    if e:
        raise ValueError('refinement zones are inconsistent:\n  - '
                         + '\n  - '.join(e))
    return True


def report():
    d = derived()
    ax = axial_plan(d)
    print('=' * 74)
    print('  Aconcagua / gmsh -- structured hex mesh plan')
    print('=' * 74)
    for k, v in [('Re_L', f"{d['Re_L']:.3e}"), ('u_tau', f"{d['utau']:.3f} m/s"),
                 ('y+ target', f'{YPLUS_TARGET:.0f}'),
                 ('first cell height y1', f"{d['y1']*1e6:.1f} um"),
                 ('BL thickness at base', f"{d['delta']*1e3:.1f} mm"),
                 ('CELLS INSIDE THE BL', f"{d['n_delta']}"),
                 ('surface cell, circumferential', f"{d['ds']*1e3:.2f} mm"),
                 ('cells around full circumference', f"{d['n_circ']}"),
                 ('butterfly cap rim', f"{d['r_cap']*1e3:.2f} mm at x = {d['x_cap']*1e3:.1f} mm"),
                 ('wall slope at the handover', f"{d['cap_angle']:.1f} deg")]:
        print(f'  {k:<34}{v}')
    print('-' * 74)
    print(f"  {'radial shell':<22}{'r_out':>9}{'cells':>7}{'ratio':>9}{'h_in':>11}{'h_out':>11}")
    for i, s in enumerate(d['shells']):
        print(f"  {'shell ' + str(i+1):<22}{s['r_out']:>9.3f}{s['n']:>7d}"
              f"{s['c']:>9.4f}{s['h_in']*1e3:>9.3f}mm{s['h_out']*1e3:>9.1f}mm")
    print(f"  {'total radial cells':<22}{'':>9}{d['n_rad']:>7d}")
    print('-' * 74)
    print(f"  {'segment':<10}{'x0':>9}{'x1':>9}{'cells':>7}{'ratio':>9}"
          f"{'h_start':>11}{'h_end':>11}")
    nx = 0
    for name, x0, x1 in segment_bounds(d):
        a = ax[name]; nx += a['n']
        print(f"  {name:<10}{x0:>9.3f}{x1:>9.3f}{a['n']:>7d}{a['c']:>9.4f}"
              f"{a['h0']*1e3:>9.2f}mm{a['h1']*1e3:>9.1f}mm")
    print(f"  {'total axial':<10}{'':>18}{nx:>7d}")
    print('-' * 74)
    core_ax = ax['up']['n'] + ax['wake1']['n'] + ax['wake2']['n'] + ax['wake3']['n']
    n_ann  = nx * d['n_rad'] * 2 * N_AZ
    n_ring = core_ax * N_RING * 2 * N_AZ
    n_core = core_ax * N_AZ * N_AZ
    print(f"  {'predicted cells: annuli':<34}{n_ann:>12,d}")
    print(f"  {'                 butterfly ring':<34}{n_ring:>12,d}")
    print(f"  {'                 butterfly core':<34}{n_core:>12,d}")
    print(f"  {'                 TOTAL (all hex)':<34}{n_ann+n_ring+n_core:>12,d}")
    print('=' * 74)
    return d, ax


if __name__ == '__main__':
    report()


# ------------------------------------------------- streamwise subdivision ---
# A transfinite quad interpolates between its boundary curves, so on a block
# whose meridian is CURVED the interior drifts off the true surface by
#     (linear blend of the end radii - true radius) x (arc - chord)
# On the LV-Haack nose that is ~1 mm across one block: three first-cell
# heights, and no amount of snapping fixes it because snapping the wall alone
# distorts the first cell.  The cure is to split the nose into blocks short
# enough that each meridian is nearly straight; the error then falls as 1/n^2.
# Stations are placed to equidistribute |r''|^(1/2), so every block carries the
# same error rather than the tip block carrying all of it.
N_NOSE_BLOCKS = 14        # -> ~7 um, about 2 % of y1
N_UP_BLOCKS   = 6         # same treatment for the upstream blend surface


def _blocks(n0, plan, xs_fn):
    """Sub-block count for a curved segment at the current H_SCALE.

    The subdivision exists to hold the interpolation error under a fraction of
    y1, and that error falls as 1/n^2 while y1 itself scales with H_SCALE, so
    n scales as 1/sqrt(H_SCALE).  The count is then reduced further if any
    sub-block would be shorter than two of its own cells.
    """
    n = max(2, int(round(n0 / math.sqrt(H_SCALE))))
    while n > 2:
        xs = xs_fn(n)
        if min(b - a for a, b in zip(xs[:-1], xs[1:])) >= \
                2.0 * max(local_h(plan, x) for x in xs):
            break
        n -= 1
    return n


def equidistribute(a, b, n, weight):
    """n+1 stations on [a, b] equidistributing `weight(x)`."""
    x = np.linspace(a, b, 20001)
    w = np.maximum(weight(x), 1e-30)
    c = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(x))])
    return np.interp(np.linspace(0.0, c[-1], n + 1), c, x)


def nose_stations(d, plan=None):
    h = 1e-6
    def w(x):
        xs = np.clip(x, h, G.L_NOSE - h)
        r2 = np.abs((G.r_body(xs + h) - 2 * G.r_body(xs) + G.r_body(xs - h)) / h ** 2)
        return np.sqrt(r2)
    f = lambda n: equidistribute(d['x_cap'], G.X_BODY_1, n, w)
    if plan is None:
        plan = axial_plan(d)['nose']
    return f(_blocks(N_NOSE_BLOCKS, plan, f))


def up_stations(d, plan=None):
    f = lambda n: np.linspace(d['x_in'], d['x_cap'], n + 1)
    if plan is None:
        plan = axial_plan(d)['up']
    return f(_blocks(N_UP_BLOCKS, plan, f))


def local_h(plan, x):
    """Local cell size of a geometric distribution at position x.

    For a geometric stack, the cell size is EXACTLY linear in distance along:
        s_k = h0 (c^k - 1)/(c - 1)  =>  h_k = h0 + (c - 1) s_k
    so subdividing a segment and re-solving each piece reproduces the parent
    distribution exactly rather than approximating it.
    """
    return plan['h0'] + (plan['c'] - 1.0) * (x - plan['x0'])


def subdivide(plan, xs):
    """Split one streamwise segment at `xs`, preserving its cell-size law."""
    out = []
    for a, b in zip(xs[:-1], xs[1:]):
        n, c = geometric(b - a, local_h(plan, a), local_h(plan, b))
        out.append(dict(x0=float(a), x1=float(b), n=n, c=c,
                        h0=local_h(plan, a), h1=local_h(plan, b)))
    return out
