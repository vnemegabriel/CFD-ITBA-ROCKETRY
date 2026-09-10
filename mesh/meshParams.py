#!/usr/bin/env python3
"""
meshParams.py -- the DEFAULT mesh parameters.  The only file you normally edit.

Every value here can also be overridden without editing, from the command
line of build.py (``--scale``, ``--sector``, ``--set KEY=VALUE``) or from a
preset in presets/.  build.py records the final parameter set next to every
.msh it writes, so a mesh is always reproducible from its own sidecar file.

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
# H_SCALE_WALL = False holds the FIRST CELL OFF EVERY WALL fixed while
# everything else scales: y1 (and therefore y+) on the lateral walls, the
# streamwise cell at the nose tip (SEGMENTS up.h_end / nose.h_start) and the
# streamwise cell off the flat base (H_WAKE_BASE).  What you want for a
# grid-convergence study with wall functions.  The tip one matters for
# quality, not just y+: the cap surface is nearly vertical at the apex and
# the core cells there are ~0.1 mm across, so a 4.5 mm streamwise cell on
# them (H_SCALE 3) pushes the boundary-face skewness past 4.
H_SCALE      = 1.0
H_SCALE_WALL = True

# H_SCALE_FIN = False does the same for the FIN: FIN_H_X (chordwise), AZ_FIN_H
# and the cell count of the two azimuthal blocks touching the fin (normal to
# it), and FIN_H_R (spanwise) keep their values while everything else scales.
# A fin is 12 mm thick with 33 / 21 mm bevels: at H_SCALE 3 with everything
# scaled it is one cell thick and its edges land two cells apart, which is
# not a fin, it is a bump.  The coarse and medium presets set this False.
H_SCALE_FIN  = True

def _h(v):
    """Scale a cell size."""
    return None if v is None else v * H_SCALE


def _hf(v):
    """Scale a FIN cell size -- unless the fin is exempt from H_SCALE."""
    return _h(v) if H_SCALE_FIN else v


def _hw(v):
    """Scale a WALL first-cell size -- unless the wall is exempt from H_SCALE."""
    return _h(v) if H_SCALE_WALL else v


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


def n_az_cells_block(k):
    """Cells in azimuthal block k (0 .. N_AZ_BLOCKS-1).

    All blocks carry n_az_cells(), except that with H_SCALE_FIN = False the
    two blocks touching the symmetry planes -- the ones the fins lie in --
    keep the unscaled N_AZ_CELLS, so the fin-normal resolution survives a
    global coarsening.  Symmetric in k <-> N-1-k, which the sector assembly
    relies on.
    """
    if not H_SCALE_FIN and k in (0, N_AZ_BLOCKS - 1):
        return max(1, int(N_AZ_CELLS))
    return n_az_cells()


def n_az_quadrant():
    """Azimuthal cells per quadrant."""
    return sum(n_az_cells_block(k) for k in range(N_AZ_BLOCKS))


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
# Radial band over which the tip closes, measured INBOARD from FIN_TIP_R.
# The method cannot make a square tip (see fin_half_thickness), so the tip is
# a chamfer and this is its height -- but only if the mesh resolves it: what
# you actually get is max(FIN_TIP_SMEAR, the radial cell at the tip), which is
# FIN_H_R when that is set.  4 mm against the 6 mm cell of `fintip` means the
# chamfer is 6 mm and the shape is whatever one cell gives.
FIN_TIP_SMEAR = 0.004        # m

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
FIN_X_LEAD   = 0.500         # m   of cylinder ahead of the root LE to include
# This one is COUPLED TO THE SWEEP, and it is the coupling that bites.  With
# FIN_EDGE_FIT the approach block runs from a plane at fin_x_start() to a
# surface that follows the leading edge, so its axial extent grows with radius
# and its cell count does not:
#
#     stretch = 1 + (FIN_TIP_LE - FIN_ROOT_LE) / FIN_X_LEAD = 1 + 0.2508 / L
#
# At the 60 mm this used to be, that is 5.1: the block is 61 mm long at the
# root and 311 mm at the tip, and its 12 uniform cells go from 5 mm to 26 mm.
# The result is a wedge of coarse cells sitting between the 5 mm cylinder and
# the 2.5 mm chord -- a 10:1 step, in the one place where the leading-edge
# shock has to be resolved.  Raising it to 500 mm brings the stretch to 1.5 and
# the cell at the tip leading edge from 25.9 mm to 7.7 mm.
#
# It is not paid for in cells, because the block is GRADED from the cylinder
# size down to FIN_H_X rather than uniform: at H_SCALE 3 the whole mesh goes
# from 530 k cells to 473 k.  What it does cost is skew, since the shear is
# now released over 500 mm instead of 61: faces above 40 deg go from 89 k to
# 116 k and the mean from 5.8 to 7.6 deg.  The MAXIMUM does not move (75.1 deg,
# at the butterfly cap) and neither does the count above 70.
#
# validate_params() computes the stretch and refuses a value that leaves it
# above FIN_LEAD_MAX_STRETCH.
FIN_LEAD_MAX_STRETCH = 2.0   # of the approach block, at the fin tip

# Spanwise: radial cell size AT THE FIN TIP.  It inserts a zone boundary AT
# FIN_TIP_R, which is what makes the fin close exactly on a node line; the
# radial stack grows from y1 to FIN_H_R over the span and continues outward
# from there.  The zone goes into the ZONE_R/ZONE_H list at the right radius --
# WAKE_ZONE_K still indexes YOUR list.  Exempt from H_SCALE when
# H_SCALE_FIN = False.
#
# None is not a neutral choice: without the node line the tip outline is
# quantised again and the fin patch comes out +4.1 % on `coarse`.
#
# It also sets the tip CHAMFER, which is max(FIN_TIP_SMEAR, the radial cell
# there) -- see fin_half_thickness.  12 mm is the working value: it matches
# `coarse` and `medium`, and the chamfer only needs to be smaller than that if
# you are after the tip vortex at incidence, which the Cd sweep is not.
# `fintip` takes it to 6 mm when you are.  Going to 1 mm costs 10.9 M annular
# cells at H_SCALE 1 against 3.0 M, and buys nothing at zero incidence.
#
# H_SCALE 6 (presets/smoke.py) cannot take a tip zone at all: scaled x6 it no
# longer fits between the tip and ZONE_R[0].  That preset declares None itself.
FIN_H_R      = 0.012          # m

# --- fin edges ON the mesh, not across it -----------------------------------
# The leading and trailing edges are swept: dx/dr = 1.56 at the LE, 0.62 at the
# TE.  A grid whose axial stations are PLANES therefore crosses them in
# diagonal, and split_symm(), which has to call a whole face wall or symmetry,
# quantises the outline to the cell.  The patch comes out as the true planform
# DILATED by one cell -- +10 % of planform area at H_SCALE 1, +21 % at 6 --
# with a serrated, zero-thickness flange ahead of the real leading edge.
#
# FIN_EDGE_FIT bends the axial stations instead: x -> x + d(x, r), a shear in
# the (x, r) plane that puts the station at FIN_LE_X_WALL exactly on x_LE(r)
# and the one at X_BODY_2 exactly on x_TE(r).  Both edges become BLOCK
# BOUNDARIES.  Nothing else changes: the mesh stays all-hexahedral and
# transfinite, and split_symm() is not touched -- its input just stops lying,
# because a node ON the leading edge has t = 0 exactly and is never moved.
#
# The price is skew, and it is not small: a face on the leading edge sits at
# 57.4 deg to the axial direction, which IS the sweep.  It cannot be reduced
# while the radial lines are circles; that would need a collar wrapped around
# the planform edge, i.e. a different topology.  Away from the two edges
# nothing goes past ~12 deg.
FIN_EDGE_FIT = True
# Radius at which the shear has died out.  None = ZONE_R[-2], the outside of
# the mid field: far enough that the blend costs ~8 deg and no more.  Pulling
# it in tightens the blend and the skew grows as 0.25 m / (r_blend - 0.2395).
FIN_EDGE_R_BLEND = None
# The trailing-edge station bulges 100 mm downstream at the tip while the base
# plane stays put, which squeezes the `tail` block to 20 % of its length there.
# This fraction of the TE shear is carried by the base station and released
# over `wake1` instead.  0 = squeeze it all into `tail`.
FIN_EDGE_TAIL_RELIEF = 0.5

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
def fin_zone_r():
    """Outer radius of the fin zone (None when FIN_H_R is off).

    FIN_TIP_R itself, so the fin closes ON a node line: the tip chamfer runs
    inboard of it and every node from there out is left where it was.
    """
    if FIN_H_R is None or not FINS_ON:
        return None
    return G.FIN_TIP_R


def shell_spec():
    """Radial shells derived from ZONE_R / ZONE_H, at call time so overrides
    of the lists are honoured.  Sizes here are UNSCALED (metres as written);
    derived() applies H_SCALE.

    Each entry: r_out, h_out, h_in (None = continue from the shell inside),
    scale ('global' or 'fin'), wake (this zone opens into the wake), name.
    The fin zone, when FIN_H_R is set, is inserted at its radius; the zones
    from ZONE_R keep their identity, so WAKE_ZONE_K indexes ZONE_R as written.
    """
    zones = [dict(r_out=r, h_out=ZONE_H[k], scale='global', wake=(k == WAKE_ZONE_K),
                  name=f'zone {k}') for k, r in enumerate(ZONE_R)]
    rf = fin_zone_r()
    if rf is not None:
        zones.append(dict(r_out=rf, h_out=FIN_H_R, scale='fin', wake=False, name='fin zone'))
        zones.sort(key=lambda z: z['r_out'])
    for k, z in enumerate(zones):
        z['h_in'] = None if k == 0 else zones[k - 1]['h_out']
        z['h_in_scale'] = None if k == 0 else zones[k - 1]['scale']
    return zones


def _hz(v, scale):
    return _hf(v) if scale == 'fin' else _h(v)

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

# ------------------------------------------------------------------ sector --
# The block structure is always built for ONE QUADRANT (theta = 0..90 deg,
# the two fins lying in its two symmetry planes).  Larger sectors are made by
# rotating copies of that quadrant and stitching them together at the
# interface planes, node for node -- see sectorAssembly.py.
#
#   'quarter'  90 deg, two symmetry planes.  Axial flow only (alpha = beta = 0).
#   'half'     180 deg, one symmetry plane (z = 0) holding two opposite fins.
#              Angle of attack in the x-y plane (alpha != 0, beta = 0).
#   'full'     360 deg, no symmetry.  Any alpha, beta, or roll; unsteady wake.
#
# Cell count scales x1, x2, x4.  Everything else -- parameters, quality,
# patch names -- is identical.
SECTOR = 'quarter'
SECTORS = {'quarter': 1, 'half': 2, 'full': 4}      # name -> number of copies

PATCHES = dict(inlet='inlet', outlet='outlet', farfield='box', symmetry='symm',
               nose='cone', body='walls', tail='tail', fins='fins')


def n_copies():
    if SECTOR not in SECTORS:
        raise ValueError(f"SECTOR must be one of {list(SECTORS)}, got {SECTOR!r}")
    return SECTORS[SECTOR]


def symmetry_fraction():
    """Fraction of the full 360 deg present in the mesh: 0.25, 0.5 or 1."""
    return n_copies() / 4.0


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
    for s in shell_spec():
        h_in = _hz(s['h_in'], s['h_in_scale']) if s['h_in'] is not None else h_prev
        h_out = _hz(s['h_out'], s['scale'])
        n, c = geometric(s['r_out'] - r_prev, h_in, h_out)
        shells.append(dict(r_out=s['r_out'], n=n, c=c, h_in=h_in, h_out=h_out,
                           wake=s['wake'], name=s['name']))
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
                ds=2.0 * np.pi * G.R_BODY / (4.0 * n_az_quadrant()),
                n_circ=4 * n_az_quadrant())


def fin_x_start():
    """Where the streamwise fin-refined block begins."""
    return G.FIN_ROOT_LE - FIN_X_LEAD


def fin_lead_stretch():
    """How much the approach block is stretched at the fin tip.

    Its upstream face is a plane at fin_x_start() and its downstream face is
    the swept leading edge, so the block is FIN_X_LEAD long at the root and
    FIN_X_LEAD + sweep long at the tip -- with the same number of cells.
    """
    if not fin_edge_fit():
        return 1.0
    return 1.0 + (G.FIN_TIP_LE - G.FIN_ROOT_LE) / FIN_X_LEAD


def fin_edge_fit():
    """Whether the fin edges are being fitted as block boundaries."""
    return bool(FIN_EDGE_FIT and FINS_ON)


def fin_edge_r_blend():
    """Radius at which the fin-edge shear has decayed to nothing."""
    if FIN_EDGE_R_BLEND is not None:
        return float(FIN_EDGE_R_BLEND)
    return float(ZONE_R[-2] if len(ZONE_R) > 1 else ZONE_R[-1])


def cyl_end():
    """Where the cylinder blocks stop and `tail` starts.

    X_BODY_2 normally.  With FIN_EDGE_FIT it is FIN_TE_X_WALL instead -- 2.5 um
    further downstream -- because that station is the one the warp bends onto
    the trailing edge, and it can only land ON the edge at every radius if it
    starts on it at the wall.  Left at X_BODY_2 the station sits 2.5 um ahead
    of its own trailing edge, every node on it keeps a half-thickness of 0.7 um
    instead of zero, and `any` then hands the whole first column of `tail`
    faces to the fin patch: 76 faces and +1.8 % of planform area, from two and
    a half microns.  Everything else is unaffected -- r_body() is evaluated at
    the actual station, so the wall stays exact, and the block either side of
    it is 2.5 um longer or shorter than it was.
    """
    return G.FIN_TE_X_WALL if fin_edge_fit() else G.X_BODY_2


def cyl_cuts():
    """Extra block boundaries inside the cylinder, as (segment name, x).

    Each one is a station the mesh will have, so each is a place the fin
    machinery can anchor to.  `cylfin` is the streamwise refinement that starts
    ahead of the root leading edge; `finchord` is the chord itself, and its
    upstream face is the one FIN_EDGE_FIT bends onto the leading edge.
    """
    cuts = []
    if FIN_H_X is not None:
        cuts.append(('cylfin', fin_x_start()))
    if fin_edge_fit():
        cuts.append(('finchord', G.FIN_LE_X_WALL))
    return sorted(cuts, key=lambda c: c[1])


def segment_bounds(d):
    """Streamwise segment end stations, in order."""
    XB = G.X_BASE
    out = [('up', d['x_in'], d['x_cap']), ('nose', d['x_cap'], G.X_BODY_1)]
    cuts = cyl_cuts()
    names = ['cyl'] + [n for n, _ in cuts]
    xs = [G.X_BODY_1] + [x for _, x in cuts] + [cyl_end()]
    out += [(names[i], xs[i], xs[i + 1]) for i in range(len(names))]
    out += [('tail',  cyl_end(), XB),
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
    out = {k: dict(h_start=_h(v['h_start']), h_end=_h(v['h_end']))
           for k, v in spec.items()}
    # first cells off the walls: the nose tip (end of `up`, start of `nose`)
    # and the flat base (start of `wake1`) follow H_SCALE_WALL, like y1
    out['up']['h_end'] = _hw(spec['up']['h_end'])
    out['nose']['h_start'] = _hw(spec['nose']['h_start'])
    out['wake1']['h_start'] = _hw(H_WAKE_BASE)
    if FIN_H_X is not None:
        hx = _hf(FIN_H_X)
        out['tail'] = dict(h_start=hx, h_end=hx)
        if fin_edge_fit():
            # GRADE the approach block from the cylinder size down to the fin
            # size, and leave `cyl` alone.  Uniform-at-FIN_H_X only made sense
            # while the block was a 60 mm strip: over the 500 mm the sweep now
            # demands it would be 100 cells of 5 mm, most of them a long way
            # from anything.  Graded it is 32, the step at either interface is
            # under 1.05 at the root, and the cells sit where the leading edge
            # is.
            out['cylfin'] = dict(h_start=out['cyl']['h_end'], h_end=hx)
        else:
            out['cyl']['h_end'] = hx
            out['cylfin'] = dict(h_start=hx, h_end=hx)
    if fin_edge_fit():
        # The chord block is the one the warp stretches and squeezes, and its
        # cell count is what sets the chordwise resolution AT EVERY SPAN
        # STATION: after the warp the axial index IS the chord fraction.
        h = _hf(FIN_H_X) if FIN_H_X is not None else out['cyl']['h_end']
        out['finchord'] = dict(h_start=h, h_end=h)
    return out


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
    tracks the spreading wake, and any zone inside it opens in proportion;
    the zones outside are cylinders.
    """
    zs = shell_spec()
    kw = wake_shell_index()
    if k > kw or x <= G.X_BASE:
        return zs[k]['r_out']
    x_out = G.X_BASE + DOWNSTREAM_L * G.L_TOTAL
    xi = min(max((x - G.X_BASE) / (x_out - G.X_BASE), 0.0), 1.0)
    rw = zs[kw]['r_out']
    grow = 1.0 + (ZONE0_R_WAKE / rw - 1.0) * xi ** WAKE_SPREAD_P
    # every zone INSIDE the wake zone (a fin zone, say) opens with it, in
    # proportion, so the layering is preserved as the band widens and no
    # inner zone gets squeezed against the blend tube at the outlet
    return zs[k]['r_out'] * grow


def wake_shell_index():
    """Index, in the derived shell list, of the zone that opens into the wake."""
    return next(k for k, z in enumerate(shell_spec()) if z['wake'])


def r_inlet():
    """Blend-tube radius at the inlet plane.  A fraction of the innermost
    zone, so it cannot outgrow the zone it lives inside."""
    return F_INLET * shell_spec()[0]['r_out']


def r_wake_out():
    """Blend-tube radius at the outlet plane: a fraction of the INNERMOST
    zone's radius there, so it cannot outgrow the zone it lives inside."""
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
    nc = n_az_cells_block(0)
    if AZ_FIN_H is not None:
        arc0 = G.R_BODY * (th[1] - th[0])
        if _hf(AZ_FIN_H) < arc0 / nc:
            k = refit(arc0, nc, _hf(AZ_FIN_H))
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
    if FIN_H_R is not None and FIN_H_R <= 0:
        e.append(f'FIN_H_R must be > 0 or None, got {FIN_H_R}')
    rf = fin_zone_r()
    if rf is not None:
        near = [r for r in ZONE_R if abs(r - rf) < 2.0 * FIN_H_R]
        if near:
            e.append(f'FIN_H_R puts a zone boundary at r = {rf:.4f} m, within two '
                     f'cells of ZONE_R entry {near[0]}: drop one of them')
        if rf >= ZONE_R[-1]:
            e.append(f'fin zone radius {rf:.4f} m is outside the farfield {ZONE_R[-1]}')
    k = WAKE_ZONE_K
    if not 0 <= k < len(ZONE_R):
        e.append(f'WAKE_ZONE_K = {k} is not a zone index')
    elif not e:
        zs = shell_spec()
        kw = wake_shell_index()
        if ZONE0_R_WAKE < zs[kw]['r_out']:
            e.append(f'ZONE0_R_WAKE {ZONE0_R_WAKE} < ZONE_R[{k}] {ZONE_R[k]}: '
                     f'zone {k} would close up downstream instead of opening out')
        if kw + 1 < len(zs) and ZONE0_R_WAKE >= zs[kw + 1]['r_out']:
            nxt = zs[kw + 1]
            e.append(f'ZONE0_R_WAKE {ZONE0_R_WAKE} >= {nxt["r_out"]} ({nxt["name"]}): '
                     f'zone {k} would swallow it at the outlet. Either lower '
                     f'ZONE0_R_WAKE or set WAKE_ZONE_K to the outermost fine zone.')
    if fin_edge_fit():
        # The warp anchors the trailing edge on the cylinder / boattail station
        # instead of giving it one of its own, which is only legitimate while
        # the two really are the same place.  They are, to 2.5 um -- but that
        # is a property of the drawing, so assert it rather than assume it.
        gap = abs(G.FIN_TE_X_WALL - G.X_BODY_2)
        if gap > 0.2 * (_hf(FIN_H_X) if FIN_H_X is not None else 0.012):
            e.append(f'FIN_EDGE_FIT anchors the fin trailing edge on X_BODY_2 = '
                     f'{G.X_BODY_2:.5f}, but the TE root is at {G.FIN_TE_X_WALL:.5f} '
                     f'({gap*1e3:.1f} mm away): the fin no longer ends on the '
                     f'boattail junction, so the TE needs a station of its own')
        cuts = [x for _, x in cyl_cuts()]
        if any(b - a < 1e-6 for a, b in zip(cuts[:-1], cuts[1:])):
            e.append(f'cylinder block boundaries collide: {cuts}. '
                     f'Lower FIN_X_LEAD or turn FIN_EDGE_FIT off')
        if not 0.0 <= FIN_EDGE_TAIL_RELIEF < 1.0:
            e.append(f'FIN_EDGE_TAIL_RELIEF must be in [0,1), got {FIN_EDGE_TAIL_RELIEF}')
        st = fin_lead_stretch()
        if st > FIN_LEAD_MAX_STRETCH:
            need = (G.FIN_TIP_LE - G.FIN_ROOT_LE) / (FIN_LEAD_MAX_STRETCH - 1.0)
            e.append(f'FIN_X_LEAD = {FIN_X_LEAD} stretches the approach block '
                     f'{st:.2f}x at the fin tip, over FIN_LEAD_MAX_STRETCH = '
                     f'{FIN_LEAD_MAX_STRETCH}: its cells there are {st:.1f} times '
                     f'the {_hf(FIN_H_X) if FIN_H_X else 0:.4g} m they are at the '
                     f'root, which puts a band of coarse cells right in front of '
                     f'the leading edge.  Use FIN_X_LEAD >= {need:.3f}, or raise '
                     f'FIN_LEAD_MAX_STRETCH to rebuild an older mesh as it was')
        rb = fin_edge_r_blend()
        if rb <= G.FIN_TIP_R + FIN_TIP_SMEAR:
            e.append(f'FIN_EDGE_R_BLEND {rb} is inside the fin tip '
                     f'{G.FIN_TIP_R + FIN_TIP_SMEAR:.4f}: nothing to blend over')
    if FIN_H_X is not None and not G.X_BODY_1 < fin_x_start() < G.FIN_ROOT_LE:
        e.append(f'FIN_X_LEAD = {FIN_X_LEAD} puts the fin block start at '
                 f'{fin_x_start():.4f}, outside the cylinder '
                 f'({G.X_BODY_1} .. {G.FIN_ROOT_LE})')
    if H_SCALE <= 0:
        e.append(f'H_SCALE must be > 0, got {H_SCALE}')
    if SECTOR not in SECTORS:
        e.append(f'SECTOR must be one of {list(SECTORS)}, got {SECTOR!r}')
    if FIN_SECTION not in ('wedge', 'diamond', 'biconvex', 'naca', 'naca_te'):
        e.append(f'unknown FIN_SECTION {FIN_SECTION!r}')
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
    if not e:
        r_prev = G.R_BODY
        for z in shell_spec():
            span = z['r_out'] - r_prev
            h = _hz(z['h_out'], z['scale'])
            if h >= span:
                e.append(f'{z["name"]}: cell size {z["h_out"]} m (scaled {h:.4g}) is not '
                         f'smaller than the zone, which is {span:.4f} m wide')
            r_prev = z['r_out']
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
        tag = s['name'] + (' (wake)' if s['wake'] else '')
        print(f"  {'shell ' + str(i+1) + '  ' + tag:<22}{s['r_out']:>9.3f}{s['n']:>7d}"
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
    est = predicted_cells(d, ax)
    print(f"  {'predicted cells: annuli':<34}{est['annuli']:>12,d}")
    print(f"  {'                 butterfly ring':<34}{est['ring']:>12,d}")
    print(f"  {'                 butterfly core':<34}{est['core']:>12,d}")
    print(f"  {'                 quadrant TOTAL':<34}{est['quadrant']:>12,d}")
    print(f"  {'sector ' + SECTOR + ' (x' + str(n_copies()) + ')':<34}{est['total']:>12,d}")
    print('=' * 74)
    return d, ax


def predicted_cells(d=None, ax=None):
    """Closed-form cell count of the quadrant, before building anything.

    Exact for the annuli and the ring; the core is an (N/2 x n_az)^2 grid per
    axial cell.  The upstream / nose sub-blocking re-solves each piece from
    the parent distribution so the axial totals are only approximate (a few
    cells per segment)."""
    if d is None:
        d = derived()
    if ax is None:
        ax = axial_plan(d)
    nx = sum(a['n'] for a in ax.values())
    core_ax = sum(ax[k]['n'] for k in ('up', 'wake1', 'wake2', 'wake3'))
    n_az_q = n_az_quadrant()                        # azimuthal cells per quadrant
    n_core_side = n_az_q // 2                       # symmetric block counts
    ann = nx * d['n_rad'] * n_az_q
    ring = core_ax * n_ring() * n_az_q
    core = core_ax * n_core_side ** 2
    q = ann + ring + core
    return dict(annuli=ann, ring=ring, core=core, quadrant=q, total=q * n_copies())


# ---------------------------------------------------------------- overrides --
# Names that build.py / presets are allowed to change.  Anything else is
# either derived or a function, and overriding it would be a silent no-op.
OVERRIDABLE = {
    'H_SCALE', 'H_SCALE_WALL', 'H_SCALE_FIN', 'FIN_H_R', 'U', 'NU', 'RHO', 'YPLUS_TARGET',
    'N_AZ_BLOCKS', 'N_AZ_CELLS', 'AZ_BLOCK_GROWTH', 'AZ_FIN_H', 'N_RING',
    'FINS_ON', 'FIN_SECTION', 'FIN_TIP_SMEAR', 'FIN_H_X', 'FIN_X_LEAD',
    'FIN_EDGE_FIT', 'FIN_EDGE_R_BLEND', 'FIN_EDGE_TAIL_RELIEF',
    'FIN_LEAD_MAX_STRETCH',
    'ZONE_R', 'ZONE_H', 'WAKE_ZONE_K', 'ZONE0_R_WAKE', 'WAKE_SPREAD_P',
    'F_INLET', 'F_WAKE_OUT', 'SEGMENTS', 'F_UP_INLET', 'F_UP_MAX',
    'H_WAKE_BASE', 'CAP_R_FRAC', 'CORE_FRAC', 'UPSTREAM_L', 'DOWNSTREAM_L',
    'X_WAKE_1', 'X_WAKE_2', 'SECTOR', 'N_NOSE_BLOCKS', 'N_UP_BLOCKS',
}


def apply_overrides(values, source='override'):
    """Set parameters from a {NAME: value} mapping, refusing unknown names.

    A typo in a preset used to be a silent no-op -- the build ran with the
    defaults and nobody noticed.  Now it is an error that names the file."""
    g = globals()
    bad = [k for k in values if k not in OVERRIDABLE]
    if bad:
        raise KeyError(f'{source}: not a mesh parameter: {bad}. '
                       f'Overridable names: {sorted(OVERRIDABLE)}')
    for k, v in values.items():
        g[k] = v


def snapshot():
    """The current parameter set, for the sidecar record written with each mesh."""
    g = globals()
    return {k: g[k] for k in sorted(OVERRIDABLE)}


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


def nose_stations(d, plan):
    h = 1e-6
    def w(x):
        xs = np.clip(x, h, G.L_NOSE - h)
        r2 = np.abs((G.r_body(xs + h) - 2 * G.r_body(xs) + G.r_body(xs - h)) / h ** 2)
        return np.sqrt(r2)
    f = lambda n: equidistribute(d['x_cap'], G.X_BODY_1, n, w)
    return f(_blocks(N_NOSE_BLOCKS, plan, f))


def up_stations(d, plan):
    f = lambda n: np.linspace(d['x_in'], d['x_cap'], n + 1)
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


if __name__ == '__main__':
    report()
