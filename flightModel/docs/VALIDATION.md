# Validation

This document gives the procedures that compare the case results with
references. Each reference has its formula and its constants, so that you
can calculate it by hand.

Do the procedures in this sequence:

1. The checks in §1, on each run.
2. Subsonic flow (§2). There are analytical references, and the
   density-based solver has the largest risk at low Mach number.
3. Transonic and supersonic flow (§3).
4. Mesh convergence (§4), before you publish a value.

## 1. Checks on each run

| Check | Correct result | Cause if not correct |
|---|---|---|
| `Cl`, `CmPitch` at α = 0, β = 0 | ≈ 0 | The mesh is not symmetric, or the run is not converged. |
| `CmRoll` at α = 0 | ≈ 0 | The fin meshes are not equal. |
| y+ on each wall patch | Most values from 20 to 100 | Change `nLayers` or `thicknessRatio`. MESH.md §4.4. |
| Minimum and maximum of T | Physically possible: near the stagnation temperature at most, not near 0 | The transient is not resolved. SOLVERS.md §8. |
| Physical time of the run | Several flow-through times (2.955 m / `Uinf`) | The run stopped before the flow was steady. SOLVERS.md §5.5. |
| `CN` from α = 2° to 6° | Linear in α | The run is not converged, or y+ on the fins is out of range. |

## 2. Subsonic flow: Barrowman

The references are the equations of Barrowman. OpenRocket uses them. The
coefficients refer to the body cross section A = πd²/4 and are per radian.

- Barrowman (1967), *The Practical Calculation of the Aerodynamic
  Characteristics of Slender Finned Vehicles*, MSc thesis, Catholic
  University of America.
- Barrowman and Barrowman (1966), *A Method for Calculating the Static
  Margin of a Slender Missile*.
- Fleeman, *Tactical Missile Design*, 2nd ed., ch. 2. Equivalent equations
  for the fins and the fin-body interference.

### 2.1 Geometry data

| Symbol | Value | Source |
|---|---|---|
| d | 0.151 m | `lRef` in `meshInfo` |
| L | 2.955 m | `lBody` in `meshInfo` |
| r | 0.0755 m | `rBody` in `meshInfo` |
| s, exposed semi-span | 0.2355 − 0.0755 = 0.160 m | `finTipR − rBody` |
| c_r, root chord | 0.301 m | [examine on the fin drawing] |
| c_t, tip chord | 0.150 m | [examine on the fin drawing] |
| l_m, length of the mid-chord line | – | [measure on the fin drawing] |
| (d_base/d)² | 0.53 | [examine on the drawing] |

> **NOTE** `finX1 − finX0` in `meshInfo` is 0.401 m. This is the x extent of
> the fin, not the root chord. Use the fin drawing for c_r, c_t and l_m.

### 2.2 Normal force and centre of pressure

| Component | `CN_α` [1/rad] | `x_cp` from the nose tip [m] |
|---|---|---|
| Nose | 2.0, for all noses of revolution | 0.466 · L_nose for an ogive |
| Cylinder | ≈ 0 in slender-body theory | – |
| Boattail | 2[(d_base/d)² − 1] = −0.94 | [calculate from the drawing] |
| Fins, N = 4, with interference | Equation below | [calculate from the drawing] |

Fin term:

```
CN_fins = K_fb · 4 N (s/d)² / (1 + √(1 + (2 l_m / (c_r + c_t))²))
K_fb    = 1 + r / (s + r)
```

### 2.3 Procedure

1. Run an angle-of-attack sweep from 0° to 6°. Refer to WORKFLOW.md §9.3.
2. Calculate the slope of `CN` against α in radians. Use the force normal
   to the body axis, not `Cl`. `Cl` is normal to the free stream.
3. Compare the slope with the sum of the components.
4. For the breakdown, add one `forceCoeffs` function object for each patch
   (`cone`, `walls`, `tail`, `fins`). Compare each term with its reference.

### 2.4 Drag at α = 0

Compare the drag components, not the total.

| Term | Reference |
|---|---|
| Skin friction | Cf · S_wet / Aref, with Cf = 0.0576 Re_L^−0.2 (turbulent flat plate, Schlichting, *Boundary-Layer Theory*, 7th ed., ch. 21) |
| Body form factor | 1 + 1/(2 · L/d), with L/d = 19.57 |
| Fin form factor | 1 + 2 t/c |
| Base drag | Hoerner, *Fluid-Dynamic Drag* (1965), ch. 3, with (d_base/d)² |

At ISA sea level:

| M | Re_L | Cf |
|---|---|---|
| 0.3 | 2.06e7 | 0.00198 |
| 0.8 | 5.50e7 | 0.00163 |
| 1.8 | 1.24e8 | 0.00139 |

The Cf correlation is for incompressible flow. Above M 0.3 it does not
include the compressibility correction.

`postProcessing/forces1` gives the pressure and viscous parts in different
columns. Compare the viscous part with the skin friction term.

## 3. Transonic and supersonic flow

No analytical reference applies to the full rocket. Use these three methods,
from the smallest to the largest work.

### 3.1 Tutorial of the template

Run the tutorial that the template comes from and get its published result:

```
$FOAM_TUTORIALS/compressible/rhoCentralFoam/biconic25-55Run35
```

The tutorial has experimental data. This method validates the solver and
the flux scheme. It separates a configuration error from a mesh error.

### 3.2 Oblique shock on a wedge

The θ-β-M relation for a wedge is exact:

```
tan θ = 2 cot β (M² sin²β − 1) / (M² (γ + cos 2β) + 2)
```

Anderson, *Modern Compressible Flow*, 3rd ed., ch. 4. The tutorial
`rhoCentralFoam/wedge15Ma5` is this case. Compare the shock angle and the
pressure ratio.

### 3.3 Cd(M) from OpenRocket

This comparison is the most important one for the trajectory. It is also
the weakest validation: OpenRocket uses semi-empirical correlations. It
finds an error of a factor of 2. It does not support a difference of 5 %.

### 3.4 Transonic drag rise

The case must show the increase of `Cd` from M 0.8 to M 1.2. For a slender
body with fins, the peak is approximately 2 to 3 times the subsonic value.
If the curve does not show it, the mesh is too coarse where the shock
occurs.

## 4. Mesh convergence

### 4.1 Procedure

1. Make three meshes. Change only `maxCellSize`, with the ratio r = √2.
   Do not change the levels. Refer to WORKFLOW.md §9.2.
2. Run the same Mach number on the three meshes.
3. Measure y+ on each mesh. The first layer thickness follows the surface
   cell size. If y+ changes much, the study contains two effects. Change
   `nLayers` to keep y+ in the same range.
4. Calculate the representative cell size h = N^(−1/3) for each mesh, with
   N the number of cells.
5. Calculate the observed order p, the extrapolated value and the grid
   convergence index (GCI) with the method of Celik et al. (2008).

### 4.2 Equations

With f₁, f₂, f₃ the results on the fine, medium and coarse mesh, and a
constant ratio r:

```
p     = ln((f₃ − f₂) / (f₂ − f₁)) / ln r
f_ext = f₁ + (f₁ − f₂) / (r^p − 1)
GCI₁₂ = 1.25 · |(f₁ − f₂)/f₁| / (r^p − 1)
```

If (f₃ − f₂)/(f₂ − f₁) is negative, the convergence is oscillatory. The
equations then do not apply.

- Celik et al. (2008), "Procedure for Estimation and Reporting of
  Uncertainty Due to Discretization in CFD Applications", *J. Fluids Eng.*
  **130**, 078001.
- Roache (1994), "Perspective: A Method for Uniform Reporting of Grid
  Refinement Studies", *J. Fluids Eng.* **116**, 405–413.

## 5. Record

For each run that you publish, keep these files with the results:

| File | Function |
|---|---|
| `<mesh-dir>/system/meshDict` | Makes the mesh again. |
| `<mesh-dir>/meshInfo` | Mesh parameters and number of cells. |
| `<mesh-dir>/log.checkMesh` | Mesh quality. |
| The STL file, or its commit | Geometry. |
| `<run-dir>/system/flowConditions` | Flight condition. |
| `<run-dir>/config.json` | Ramp and stop criterion. |
| y+ for each patch from `log.rhoCentralFoam` | Wall resolution. |

With these files, another person can make the same result again.
