# Aconcagua — gmsh structured mesh: design notes

*Why the mesh is built the way it is.  For how to use it see
[WORKFLOW.md](WORKFLOW.md) (Spanish); every parameter is listed in
[PARAMETERS.md](PARAMETERS.md); half and full meshes in
[SECTORS_AND_AOA.md](SECTORS_AND_AOA.md).*

3,538,188 cells per quadrant, 100 % hexahedra, fins resolved, amplified wake.
checkMesh: non-orthogonality max 86.7 deg, mean 5.35, 7,563 faces above 70 of
10.7 million; skewness max 2.75.  Half (2 quadrants) and full (4) meshes are
assembled from this quadrant with identical quality.

---

## 1. Mesh zones

**Streamwise segments and radial levels.** Both axes are piecewise-linear so every
block is visible — the picture is not to scale.

![meridional zones](img/zones_meridional.png)

The copper line is **shell 1**, the fine radial band. Downstream of the base it
**opens out from 0.35 m to 0.70 m** so it tracks the spreading wake instead of the
wake growing out of it. Plum is the butterfly core and ring, which exist only where
the axis is fluid — upstream of the nose cap and downstream of the base.

**Cross-plane tiling** of the quarter: 12 azimuthal blocks × 3 cells, block widths
graded toward both symmetry planes because that is where the fins are. The butterfly
core is a 6 × 6 grid of sub-blocks: a single quad core can present only two edges to
the ring, so `N_AZ_BLOCKS > 2` requires an (N/2) × (N/2) core.

![cross-plane zones](img/zones_crossplane.png)

| Level | Zone | Extent | Cells | Present where |
|---|---|---|---|---|
| L0 | core | butterfly square, half-width 0.45·R_i | 18 × 18 | axis fluid only |
| L1 | ring | square → inner circle R_i(x) | 10 radial | axis fluid only |
| L2 | zone 0 | wall → `ZONE_R[0]` = 0.35 m (0.70 m at the outlet) | 57 radial | everywhere |
| L3 | zone 1 | → `ZONE_R[1]` = 2.00 m | 21 radial | everywhere |
| L4 | zone 2 | → `ZONE_R[2]` = 11.775 m — **this is the farfield** | 14 radial | everywhere |

| Segment | x₀ | x₁ | Axial cells | Wall patch |
|---|---|---|---|---|
| up (×6) | −17.730 | 0.026 | 79 | — |
| nose (×14) | 0.026 | 0.800 | 153 | cone |
| cyl | 0.800 | 2.469 | 209 | walls |
| cylfin | 2.469 | 2.830 | 72 | walls |
| tail | 2.830 | 2.955 | 25 | tail |
| wake1 | 2.955 | 3.255 | 110 | — |
| wake2 | 3.255 | 6.955 | 222 | — |
| wake3 | 6.955 | 41.370 | 91 | — |

**Fins.** Planform and section:

![fin planform and section](img/zones_fin.png)

---

## 2. Parameters — `meshParams.py`

You give **cell sizes**; counts are derived from `c = (L−h₁)/(L−h_N)` and
`n = 1 + ln(h_N/h₁)/ln c`. Both closed form.

`meshParams.py` holds the defaults; anything in it can be overridden from
`build.py` (`--preset`, `--scale`, `--sector`, `--set KEY=VALUE`) without
editing, and the parameter set actually used is written next to every mesh as
`<name>.params.py`.  The tables below describe the defaults.

### Global coarsity — one number

```python
H_SCALE      = 1.0     # multiplies every cell size, divides every cell count
H_SCALE_WALL = True    # False holds y1 (and y+) fixed while the rest scales
```

Everything re-solves from the scaled sizes, so distributions stay correct rather
than being thinned out. `build.py --scale <H>` sets it from the command line (the legacy `AG_COARSE=<n>` environment variable still works).

| `H_SCALE` | cells | non-orth mean | max skew | build |
|---|---|---|---|---|
| 0.8 | 7,430,784 | — | — | |
| **1.0** | **3,538,188** | **1.49°** | 3.16 | ~2 min |
| 2.0 | 590,304 | 1.71° | 2.72 | ~25 s |
| 3.0 | 124,716 | 2.02° | — | ~10 s |
| 4.0 | 70,596 | 2.26° | — | |
| 6.0 | 30,768 | 2.73° | — | ~5 s |

Set `H_SCALE_WALL = False` for a grid-convergence study: y⁺ stays at 32 while every
other direction coarsens.

### Refinement zones — the block to edit

```python
ZONE_R = [0.35, 2.00, 11.775]    # where each zone ENDS; LAST ONE IS THE FARFIELD
ZONE_H = [0.020, 0.200, 1.600]   # cell size at the OUTER edge of each zone
WAKE_ZONE_K   = 0                # which zone opens into the wake
ZONE0_R_WAKE  = 0.70             # that zone's outer radius at the outlet
WAKE_SPREAD_P = 0.5              # turbulent wake spreads as x^(1/2)
F_INLET    = 0.63                # blend tube at the inlet, as a FRACTION of zone 0
F_WAKE_OUT = 0.50                # blend tube at the outlet, same
```

`R_FAR` and `SHELLS` are **derived** from `ZONE_R`/`ZONE_H` — the farfield radius is
stated once, as `ZONE_R[-1]`. The blend tubes are fractions rather than metres so they
cannot outgrow the zone that contains them.

`validate_params()` runs on every build and names what is wrong:

| Rule | Error if broken |
|---|---|
| `ZONE_R` strictly increasing | `ZONE_R must increase: [...]` |
| `ZONE_R[0] > R_BODY` | `ZONE_R[0] = ... is inside the body` |
| `ZONE_R[k] ≤ ZONE0_R_WAKE < ZONE_R[k+1]`, k = `WAKE_ZONE_K` | that zone would close up, or swallow the next |
| `FIN_X_LEAD` keeps the fin block start on the cylinder | names the resulting x and the valid range |
| `H_SCALE > 0` | out of range |
| `ZONE_H[k]` smaller than zone k | `ZONE_H[k] is not smaller than zone k, which is ... wide` |
| `0 < F_INLET, F_WAKE_OUT < 1` | out of range |
| `CORE_FRAC < 1/√2` | the core square would poke through its own ring |
| `N_AZ_BLOCKS` even and ≥ 2 | the core grid is undefined otherwise |

### Everything else

| Parameter | Value | What it does |
|---|---|---|
| `N_AZ_BLOCKS` | 12 | azimuthal blocks per quadrant, **must be even** (the core is an (N/2)² grid) |
| `N_AZ_CELLS` | 3 | cells per azimuthal block → 36 per quadrant, 144 around |
| `AZ_BLOCK_GROWTH` | 1.50 | block-width ratio, from the symmetry planes inward |
| `AZ_FIN_H` | 5e-4 m | first azimuthal cell at r = R_BODY in the block touching each symmetry plane. `None` = uniform |
| `N_RING` | 10 | cells across the butterfly ring |
| `YPLUS_TARGET` | 32 | y⁺ at the first cell **centre**; y₁ = 303 µm, 31 cells inside δ |
| `SEGMENTS` | `h_start, h_end` | streamwise cell size. `None` = continue from the neighbour |
| `F_UP_INLET` | 0.50 | inlet streamwise cell, as a fraction of the first upstream sub-block (capped by `F_UP_MAX`) |
| `H_WAKE_BASE` | 5e-4 m | first cell off the flat base — it is a wall |
| `CAP_R_FRAC` | 0.10 | butterfly cap rim / R_BODY. Rim 7.55 mm at x = 26.3 mm |
| `CORE_FRAC` | 0.45 | core half-width / R_i. Must stay < 0.707 (validated) |
| `X_WAKE_1` / `X_WAKE_2` | 0.30 / 4.00 m | near wake runs to 27 body diameters |
| `WAKE_ZONE_K` | 0 | which zone opens into the wake |
| `FIN_H_X` / `FIN_X_LEAD` | 0.005 / 0.060 m | streamwise refinement over the fin chord |
| `H_SCALE` / `H_SCALE_WALL` | 1.0 / `True` | global coarsity |
| `UPSTREAM_L` / `DOWNSTREAM_L` | 6 / 13 | domain, in body lengths |
| `FINS_ON` | `True` | |
| `FIN_SECTION` | `'wedge'` | `wedge` / `diamond` / `biconvex` / `naca` / `naca_te` |
| `FIN_TIP_SMEAR` | 0.004 m | radial band closing the tip taper |

Geometry constants live in `aconcaguaGeom.py`: five numbers define the body
(`L_NOSE`, `L_CYL`, `L_TAIL`, `R_BODY`, `R_BASE`), everything else derived, and
`validate()` asserts the invariants.

## 3. Scripts

```
python build.py --plan                    # print the plan before building anything
python build.py --preset fine             # build -> snap -> fins -> audit -> write   (~2 min)
python build.py --preset smoke            # pipeline test, ~25 s
python build.py --preset medium --sector full
cd <run>; ./Allmesh <file.msh>            # gmshToFoam -> meshInfo -> patch types -> checkMesh -> renumberMesh
```

| File | What it is |
|---|---|
| `aconcaguaGeom.py` | geometry of record, replaces the STL. Self-verifying. |
| `meshParams.py` | defaults, validation, the cell-count solver.  The only file you normally edit |
| `presets/*.py` | named parameter sets layered on the defaults |
| `build.py` | the command line: presets, overrides, sector, output, sidecars |
| `blockTools.py` | memoising structured-block layer over gmsh |
| `buildHexBody.py` | the block topology |
| `meshIO.py` | gmsh model -> numpy arrays; .msh v2.2 writer; `.meshInfo` / `.params.py` sidecars |
| `meshFinish.py` | the quadrant pipeline: build, mesh, extract, snap, deform, classify |
| `finPatch.py` | fin deformation and symmetry/fin face classification, on arrays |
| `sectorAssembly.py` | quadrant -> half / full by rotate-and-stitch |
| `meshQuality.py` | OpenFOAM's quality measures, computed before OpenFOAM does |
| `../case/fixPatchTypes.py` | **run right after gmshToFoam, not optional**; types from `meshInfo` |

Everything after `gmsh.model.mesh.generate(3)` works on plain arrays: no
per-node gmsh API calls, one vectorised pass each for the snap and the fin
deformation, and a writer that emits exactly what gmshToFoam consumes.  Only
nodes referenced by a hexahedron are written, so spline control points never
reach the file.

Setup:

```
pip install gmsh numpy scipy
sudo apt install libglu1-mesa     # Linux/WSL: import gmsh fails with OSError: libGLU.so.1 without it
```

The quadrant `.msh` is 520 MB ASCII — **generate it locally, do not copy it**;
the `.params.py` sidecar rebuilds it.

`gmshToFoam` makes every patch `type patch`. `symm` fails loudly; **`cone`/`walls`/
`tail`/`fins` fail silently** — wall functions, yPlus and forceCoeffs all quietly
wrong. `fixPatchTypes.py` sets them from the `patchTypes` entry of `meshInfo` and
errors on anything unrecognised.

Patches: `inlet` `outlet` `symm` `box` `cone` `walls` `tail` `fins` — 266,244 faces
in the quadrant, exactly the boundary-face count (build.py checks this and refuses
to write otherwise), so nothing lands in `defaultFaces`.

---

## 4. Additional modifications

### How the fins went in

No blocks were added. The azimuthal coordinate is **deformed**:

```
theta -> theta_f(x,r) + theta * (90 - 2 theta_f) / 90     theta_f = arcsin(t_half / r)
```

A node that was on the symmetry plane lands at `z = -t_half`, which **is** the fin
surface, exact to **0.0 nm**; wetted area 718.62 cm² against 723.00 cm² analytic.

It works only because the real fin is bevelled — thickness goes to zero continuously
at the leading and trailing edges, so the deformation relaxes to nothing there. Both
symmetry planes get it: the quarter contains two half fins, one lying in each.
`theta_f` is 2.28° at the root, 0.73° at the tip.

**Which faces are fin.** A symmetry-plane face with *any* displaced node is off the
plane, so it is tagged `fins`; only faces whose four nodes stayed put remain `symm`.
The patch therefore carries a one-cell rim around the planform, which shrinks
first-order with `FIN_H_X` and `FIN_H_R`: against the exact wetted area of the
deformed surface, 376.0 cm² per half-fin, the mesh gives +3.8 % at 5/12 mm,
+2.0 % at 2.5/6 mm and +0.9 % at 1.25/3 mm.  The **normal** is right everywhere —
the fin is a flat plate in the plane — so this is not the staircase that corrupts
wall shear stress; it is a quantised outline.  Note that the exact reference is
376.0 and not the 361.5 cm² nominal planform: the nominal figure omits the tip
smear band and treats the bevels as flat, and comparing against it inflated the
apparent discretisation error by about 4 %.  The same test decides, when two
quadrants are stitched into a half or full mesh, which interface faces become
interior and which stay as the two walls of the now full-thickness fin — one rule,
so the quarter, half and full patches cannot disagree (see `sectorAssembly.py`).

### Refining around the fins

Three directions, three places to edit.

| Direction | Parameter | Now | Effect |
|---|---|---|---|
| **azimuthal** | `AZ_FIN_H` | 5e-4 m | first cell off the fin surface → y⁺ ≈ 53 |
| | `AZ_BLOCK_GROWTH` | 1.50 | how fast blocks widen away from the fin |
| **streamwise** | `FIN_H_X` | 0.005 m | cell over the fin chord; splits `cyl` into `cyl` + `cylfin` and refines `tail`. `None` = off |
| | `FIN_X_LEAD` | 0.060 m | cylinder included ahead of the root leading edge |
| **radial** | `ZONE_R` / `ZONE_H` | see below | add a zone just outside the tip |

What those settings currently produce:

| | value |
|---|---|
| first cell off the fin, azimuthal | 0.500 mm |
| streamwise over the chord, `cylfin` x 2.4692 … 2.8300 | 5.000 mm, 72 cells |
| streamwise over the boattail, `tail` | 5.000 mm, 25 cells |
| radial at the fin tip, r − R = 0.16 m | 11.89 mm |
| fin faces | 5,218 |

**Radial refinement at the tip** is the one that needs a zone rather than a stack
parameter. The tip sits at r = 0.2355 m, inside zone 0, where the radial stack has
already grown to ~12 mm. `FIN_H_R` inserts a zone boundary at the outer edge of the
tip smear (`FIN_TIP_R + FIN_TIP_SMEAR` = 0.2395 m, so the fin closes on a node line)
with that radial size:

```python
FIN_H_R = 0.006          # presets/fintip.py
```

That gives 85 cells from the wall to the tip at ratio 1.036, 6 mm at the tip instead
of 12, and 129 radial cells in total against 92. The zone is inserted into the
ZONE_R list at its radius, `WAKE_ZONE_K` keeps indexing the user's list, and every
zone inside the wake zone opens with it in proportion downstream of the base — an
inner zone that stayed cylindrical while the wake tube grew to 0.35 m at the outlet
was squeezed to negative width, 816 inverted cells, which `build.py` refused to
write. A transfinite block carries its radial count everywhere, so the fin zone's
cells also run the length of the domain: that is the price of a fin zone, and why
`fine` does not have one by default.

**Coarsening and the fin.** `H_SCALE_FIN = False` (the coarse and medium presets)
exempts the fin from `H_SCALE`: `FIN_H_X`, `AZ_FIN_H`, the cell count of the two
azimuthal blocks touching the fin, and `FIN_H_R`. With everything scaled by 3 the
12 mm fin was one cell thick with its edges two cells apart — a bump with a fin's
planform — and the `fins` patch carried a 30 % rim. The per-block azimuthal count
is consistent because a core u-edge shares its count with the ring block it faces
(`azu{i}` ↔ block N−1−i) and symmetric in k ↔ N−1−k, which the sector stitching
needs.

That clustering reaches the nose apex through the butterfly core (the core's
edges share the ring blocks' azimuthal distributions), so the cells on the axis at
the tip are ~0.1 mm across.  With the streamwise cell at the cap scaled to 4.5 mm
the cap's boundary faces there reached skewness 4.6 (checkMesh flags > 4).
`H_SCALE_WALL = False` therefore also holds the streamwise first cell at the tip and
off the base, not just `y1`: 2.5 in the coarse preset against 2.75 in `fine`.

**Audit vs checkMesh.** `meshQuality.py` now uses OpenFOAM's own skewness
definitions (internal faces normalised by the face extent in the skew direction,
boundary faces against the owner's normal projection); an earlier version
normalised by √area and over-reported stretched faces by 2–3×.  Non-orthogonality
still differs (vertex-averaged vs volume-weighted cell centres): 75.9° here against
86.7° in checkMesh on the same mesh.  checkMesh's numbers are the ones that count.

### Azimuthal quality

| non-orth mean | **1.49°** |
|---|---|
| non-orth max | 75.85° |
| faces > 70° | 6,999 of 10,481,442 |
| max skewness | 3.16 |
| mean skewness | 0.012 |

Two things produce this, and both are needed.

**Radial curves are canonicalised inward → outward before creation.** A radial curve is
edge 0 of one annulus block and edge 2 of its neighbour, so the two blocks request it in
opposite directions, and whichever call reaches the memo first fixes the direction gmsh
runs the progression in. Without canonicalisation one azimuthal column gets the intended
wall-clustered stack and every other column gets its reciprocal — first radial cell
11.5 mm against 30.5 mm at the same station, cells shrinking outward instead of growing.
The same rule already applies to the azimuthal families in `_az`.

**Small angular spans.** A transfinite quad blends radius linearly in index, so the
interpolation error scales with the block's arc-to-chord sagitta. The widest sector is
now 16.4° (sagitta 0.0103) against 45° (0.0761).

### Still open

checkMesh counts 7,563 faces above 70° out of 10.7 million (the vertex-centred audit
in `meshQuality.py` says 6,999 and a 75.9° maximum; checkMesh's volume-weighted
centroids give 86.7°).  With `-allGeometry` it also flags 3,030 cells of aspect
ratio up to 3,026: all in the farfield (r > 2 m) just behind the base, where the
0.5 mm first streamwise cell off the base meets 1.6 m radial cells — a transfinite
block carries one axial distribution across all radii.  Harmless there; removing it
would need a radially varying axial distribution, i.e. a different topology.
Neither is blocking; `fvSchemes` carries `limited corrected 0.33` and
`nNonOrthogonalCorrectors 1`.

### Toolchain constraints the code depends on

Each of these forces something in `blockTools.py` or `buildHexBody.py`. Changing that
code without honouring them produces a mesh that builds and is wrong.

| Constraint | What it forces |
|---|---|
| a memoised shared entity is created once, in whichever direction the first caller asks | canonicalise direction before creating any curve — radial, azimuthal, axial |
| gmsh spaces transfinite points on a spline by **arc length, not parameter** | control points cannot dictate node positions; sample densely for curve accuracy and pass the real progression |
| a power blend `ξ^q` with 1 < q < 2 has **unbounded curvature** at ξ = 0 | the upstream blend is tangent-matched, `Rᵢ² = r_cap² + 2 r_cap s_cap t + λt²` |
| a transfinite quad drifts off a **curved** meridian by (linear blend of end radii − true radius) × (arc − chord) | the nose is 14 blocks with stations equidistributing \|r″\|^½; error falls as 1/n² |
| every gmsh model Point carries a mesh node | spline control points are un-meshed by `drop_control_nodes()` |
| a gmsh physical group is **per-surface**, but one block face carries both fin and symm | the fin patch is split at element level in the `.msh` |
| `c ** n` overflows once n·ln c passes ~709 | `stack_sum` is log-guarded |
| `gmshToFoam` types **every** patch `patch` | `fixPatchTypes.py` must list every patch and skip the `FoamFile` header |

`Aref` references **the body cross-section**, not the fin semi-span, times the
fraction of 360° in the mesh: `0.0044770 m²` for the quadrant, `lRef = 0.1510 m`.
Both, and the wall patch list, are written by `build.py` into the `.meshInfo`
sidecar that `controlDict` includes — `forceCoeffs` aborts on a patch name it
cannot find, and a hand-maintained list drifted from the mesh more than once.

## Current quality

checkMesh -allGeometry -allTopology on the fine quadrant (OpenFOAM v2412), against
the snappyHexMesh mesh it replaced:

| | snappyHexMesh | this mesh (quadrant) |
|---|---|---|
| cells | 2,475,822 | 3,538,188 |
| hexahedra | 96.9 % | **100 %** |
| polyhedra / prisms | 67,047 / 10,045 | **0 / 0** |
| concave cells | 6,026 | **0** |
| illegal faces | 11 | **0** |
| cells with volume ≤ 0 | — | **0** |
| cells inside δ | 12–20 | **31** |
| wall residual | — | **0.00 nm** |
| max skewness | 2.51 | 2.75 |
| max aspect ratio | 24.7 | 3,026 (3,030 farfield cells; 63 elsewhere) |
| non-orth mean | 6.81° | **5.35°** |
| non-orth max | — | 86.7° |
| non-orth > 70° | 496 | 7,563 |
| checkMesh | — | passes the default checks; fails 3 of `-allGeometry` (aspect ratio, determinant, interpolation weight) |
| y⁺ on the fins | — | ~53 |

The half and full meshes reproduce these numbers exactly, per copy.
