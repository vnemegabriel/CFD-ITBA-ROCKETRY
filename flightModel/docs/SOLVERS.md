# Solver and models

This document gives the dictionaries of the case templates `case-central/`
and `case-pimple/`, the models that they select, and the references. Use it to calculate the
case values by hand and to compare them with the case.

> **CAUTION** The template is not validated. Do not use the results before
> you do the procedures in [VALIDATION.md](VALIDATION.md).

## 1. Summary

| Item | Value |
|---|---|
| Templates | `case-central/` (M > 1, started from pimple, §17), `case-pimple/` (M < 1, §16) |
| Solvers | `rhoCentralFoam`; `rhoPimpleFoam` |
| Formulation | Compressible. Central: density-based, explicit. Pimple: pressure-based, `transonic` |
| Time | `localEuler` in both, time in iterations (§5). Central: Courant ramp (§5.2) |
| Fields in `0.orig/` | `U`, `p`, `T`, `k`, `omega`, `nut`, `alphat` |
| Turbulence | k-ω SST with wall functions |
| Gas | Air, perfect gas, constant Cp, Sutherland viscosity |
| Input | `Minf`, `pInf`, `Tinf`, `alpha`, `beta`, `Ti`, `nuRatio` |
| Initial condition | Uniform free stream |

Two solvers, split at M = 1. The central-upwind flux of `rhoCentralFoam`
scales its numerical dissipation with the speed of sound, not with the flow
speed, so it is too dissipative at low Mach. The OpenFOAM ToolChain report
(Martí Valverde and Stemmer, 2026, doi:10.5281/zenodo.19823943, §6) measured
a CA error of about 100 % against wind-tunnel data with `rhoCentralFoam` at
low Mach, and better results with `rhoPimpleFoam` for M < 0.8. Above
M 1.5 `rhoPimpleFoam` gave nonphysical results at the nose tip.

The price is a change of discretization at M = 1: a step in the Cd–Mach
curve there is not separable from the physics. To measure it, run one Mach
number with both solvers on the same mesh (§16.6).

The models, the gas, the boundary conditions, the function objects and the
flight condition are the same files in both templates.

## 2. Files of the template

This section lists `case-central/`. `case-pimple/` holds only the files that
differ; §16.1 lists them.

| File | Edit | Function |
|---|---|---|
| `system/flowConditions` | yes | Flight condition. §3. |
| `config.json` | yes | Parameters of the Courant ramp and of the stop criterion. §5. |
| `system/flowDerived` | no | Values calculated from `flowConditions`. §3. |
| `0.orig/*` | no | Initial and boundary conditions. §9. |
| `constant/thermophysicalProperties` | no | Gas model. §7. |
| `constant/turbulenceProperties` | no | Turbulence model. §6. |
| `constant/meshInfo` | no | `newCase.sh` replaces it with the `meshInfo` of the mesh. MESH.md §7. |
| `system/controlDictBase.j2` | yes, with care | Time control. §5, §10. |
| `system/FOs/FOcommon` | yes, with care | Function objects of both templates, except `residuals`. §10. |
| `system/controlDict` | no | Includes `controlDictBase` and adds `maxCo 0.01`. |
| `system/fileUpdaterControlDict/controlDict_{1..7,final}.j2` | no | Steps of the Courant ramp. §5. |
| `system/fileUpdater.j2` | no | Schedule of the Courant ramp. §5. |
| `system/FOs/FOrunTimeControl.j2` | no | Stop criterion. §5. |
| `system/fvSchemes`, `fvSchemesUpwind` | no | First-order schemes for the start. §11. |
| `system/fvSchemesVanLeer` | yes, with care | Second-order schemes. §11. |
| `system/fvSolution` | no | Linear solvers. §12. |
| `system/decomposeParDict` | through `--np` | Parallel decomposition. §13. |
| `system/topoSetDict`, `system/refineMeshDict` | `margin`, `wake`, `tipBand` only | Local refinement. WORKFLOW.md §7. |
| `Allclean` | no | Removes the mesh, the fields and the logs. |

> **CAUTION** Use `foamDictionary -set` only on `system/flowConditions`.
> On a file with `#include` or `#eval`, `foamDictionary -set` writes the
> expanded values. The file then does not follow `flowConditions`.

## 3. Flight condition

### 3.1 flowConditions

| Entry | Default | Unit | Function |
|---|---|---|---|
| `Minf` | 1.80 | – | Free-stream Mach number. |
| `alpha` | 0 | deg | Angle of attack. The velocity turns from +x to +y. |
| `beta` | 0 | deg | Angle of sideslip. The velocity turns to +z. |
| `pInf` | 101325 | Pa | Static pressure, ISA sea level. |
| `Tinf` | 288.15 | K | Static temperature, ISA sea level. |
| `Ti` | 0.005 | – | Turbulence intensity. |
| `nuRatio` | 5 | – | Ratio `nut/nu` in the free stream. |

To change a value, use `newCase.sh` (WORKFLOW.md §4) or this command:

```bash
foamDictionary system/flowConditions -entry Minf -set 1.5
```

### 3.2 flowDerived

`flowDerived` includes `flowConditions` and calculates the values below with
`#eval`. The fields in `0.orig/` and `controlDictBase` include `flowDerived`.
Thus, a change to `flowConditions` changes all the conditions.

| Value | Formula | M 0.3 | M 0.8 | M 1.8 |
|---|---|---|---|---|
| `aInf` [m/s] | √(γ R T) | 340.29 | 340.29 | 340.29 |
| `Uinf` [m/s] | M · a | 102.09 | 272.23 | 612.53 |
| `rhoInf` [kg/m³] | p / (R T) | 1.225 | 1.225 | 1.225 |
| `muInf` [Pa s] | As √T / (1 + Ts/T) | 1.790e-5 | 1.790e-5 | 1.790e-5 |
| `nuInf` [m²/s] | μ / ρ | 1.461e-5 | 1.461e-5 | 1.461e-5 |
| `kInlet` [m²/s²] | 1.5 (U · Ti)² | 0.391 | 2.779 | 14.07 |
| `nutInlet` [m²/s] | nuRatio · ν | 7.307e-5 | 7.307e-5 | 7.307e-5 |
| `omegaInlet` [1/s] | k / νt | 5349 | 3.803e4 | 1.925e5 |
| `alphatInlet` [kg/(m s)] | ρ νt / Prt | 1.053e-4 | 1.053e-4 | 1.053e-4 |

The values are for ISA sea level, `Ti 0.005`, `nuRatio 5`.

| Value | Formula |
|---|---|
| `Ux`, `Uy`, `Uz` | U·(cos α cos β, sin α cos β, sin β) |
| Drag direction `dx`, `dy`, `dz` | (cos α cos β, sin α cos β, sin β) |
| Lift direction `lx`, `ly`, `lz` | (−sin α, cos α, 0) |

Constants in `flowDerived`: `Rgas 287.05`, `gammaGas 1.4`, `Prt 0.85`. They
must agree with `thermophysicalProperties` (§7). The Sutherland constants
are written in the `muInf` formula.

## 4. Equations and numerical method

`rhoCentralFoam` does not solve a pressure equation. It advances the
conserved variables ρ, ρU and ρE with a central-upwind flux. Then it applies
the diffusive terms implicitly. Thus `fvSolution` has no `SIMPLE` or
`PIMPLE` dictionary and no relaxation factors.

**Flux.** `fluxScheme Kurganov`. The scheme divides the face flux with the
local wave speeds a⁺ and a⁻. It captures shocks without a Riemann solver and
without artificial viscosity.

- Kurganov and Tadmor (2000), "New High-Resolution Central Schemes for
  Nonlinear Conservation Laws and Convection-Diffusion Equations",
  *J. Comput. Phys.* **160**, 241–282.
- Kurganov, Noelle and Petrova (2001), *SIAM J. Sci. Comput.* **23**(3),
  707–740.
- Greenshields, Weller, Gasparini and Reese (2010), "Implicit and explicit
  schemes for flows of aerodynamic and turbomachinery fluids", *Int. J.
  Numer. Meth. Fluids* **63**, 1–21. Implementation and validation in
  OpenFOAM.

**Reconstruction.** The primitive variables go to the faces with the van
Leer limiter: `reconstruct(rho) vanLeer`, `reconstruct(U) vanLeerV`,
`reconstruct(T) vanLeer`. At a discontinuity the limiter decreases the
scheme to first order. In smooth regions the scheme stays second order.
`vanLeerV` limits in the direction of the largest gradient. van Leer (1979),
*J. Comput. Phys.* **32**, 101–136.

## 5. Local time step, Courant ramp and stop criterion

### 5.1 Local time step

Both templates use `ddt localEuler`, as the OpenFOAM ToolChain does (report
§5.3). Each cell advances with its own time step, the step that gives
`maxCo` in that cell. The run is a pseudo-transient march to the steady
state: the intermediate states are not physical, the converged state is.

`deltaT 1` in `controlDict`, so the time counts iterations: `endTime`,
`writeInterval`, the ramp and the stop criterion are all in iterations.

On `h1` at M 1.8 the local step goes from 1.5e-10 s in the smallest cells to
1.6e-5 s in the far field. With the global `Euler` step, which the template
used before, all the domain advanced at 1.5e-10 s: about 6.8e5 steps per
flow-through time, about two weeks per flow-through time on 8 cores.

| Entry | Central: read from | Pimple: read from | Value |
|---|---|---|---|
| `maxCo` | `controlDict` (the ramp swaps it) | `fvSolution` `PIMPLE` | ramp to `final_Co`; pimple `maxCo` in `config.json` |
| `maxDeltaT` | `controlDict` | `fvSolution` `PIMPLE` | 1 s, no limit |
| `rDeltaTSmoothingCoeff` | `controlDict` | `fvSolution` `PIMPLE` | 0.1 |
| `rDeltaTDampingCoeff` | not read | `fvSolution` `PIMPLE` | pimple 0.9 |

Source: `applications/solvers/compressible/rhoCentralFoam/setRDeltaT.H` and
`rhoPimpleFoam/setRDeltaT.H` of v2412.

> **CAUTION** With `localEuler`, the "Mean and max Courant Numbers" line of
> `rhoCentralFoam` uses `deltaT` of `controlDict` (1), not the local step.
> Values of order 1e7 are normal. `rhoPimpleFoam` does not write the
> Courant number. Look at `Flow time scale min/max` instead: the minimum is
> the step of the most restrictive cell. A sudden drop shows a cell that
> diverges.

### 5.2 Courant ramp (central only)

A uniform free stream on a no-slip wall causes a large transient. The ramp
of the ToolChain (report §5.4.1, `templates/rhoCentralFoam`) starts with
`maxCo` 0.01 and first-order reconstruction, and doubles `maxCo` step by
step. The `timeActivatedFileUpdate` function objects in `system/fileUpdater`
copy a file on `system/controlDict` or on `system/fvSchemes` at each switch
iteration. `runTimeModifiable yes` makes the solver read the new file.

| Phase | Iterations | `maxCo` |
|---|---|---|
| Upwind | 0, 20, 35, 50, 65, 80, 100 | 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.5 |
| van Leer | `vanLeer` + 0, 20, 35, 150, 300, 800, 2000 | 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.5 |
| Final | `minIter`, or the end of the ramp + 1 if later | `final_Co`, with `runTimeControl` |

Each `maxCo` is limited to `final_Co`. The final step never comes before the
end of the ramp: `rhoCentralFoam` has no damping of the local step, and a
jump of `maxCo` from 0.04 to 0.5 gave a floating point exception in the
thermophysical model six iterations later (`h1`, M 1.8). With the default
values the ramp ends at iteration 4000, so the final step is at 4001.

`rhoPimpleFoam` has no ramp: the ToolChain ran it without one, and it was
stable from the uniform free stream on `h1` at M 0.3 (§16.6).

### 5.3 config.json

`newCase.sh` renders all `*.j2` files with these values.

| Key | Central | Pimple | Function |
|---|---|---|---|
| `maxCo` | – | 0.5 | Courant number of each cell. |
| `vanLeer` | 2000 | – | Iteration of the change from upwind to van Leer. `-1` keeps upwind. |
| `minIter` | 3000 | 1000 | Earliest iteration of the stop criterion. Central: change to `controlDict_final`, not before the end of the ramp (§5.2). |
| `final_Co` | 0.5 | – | Maximum `maxCo` of the ramp. |
| `endIter` | 20000 | 20000 | Last iteration if the criterion does not stop the run. |
| `coeffs_fields` | `["Cd"]` | `["Cd"]` | Coefficients in the stop criterion. |
| `coeffs_variation` | `[0.007]` | `[0.007]` | Absolute tolerance of each coefficient: 0.007 is 2 % of `Cd` = 0.35, that of `m08-a00-h0`. A higher `Cd` gets less than 2 %, a lower one more. |
| `coeffs_window` | `[1000]` | `[1000]` | Window of each coefficient, in iterations. |

The sweep stops on `Cd` alone, also at α = 4 and 8°: check that `Cl` and
`CmPitch` are flat in `coefficient.dat` too. A tolerance for `CmPitch` has
not been chosen: about the nose tip with `lRef` = 0.151 m it is of order
10 to 100 times `Cl`. You can add `Cl` and `CmPitch` to `coeffs_fields`. The tolerance
is absolute, so it also applies to a coefficient with an average near 0.

The values come from the ToolChain, which used another rocket and other
meshes. Adjust them from the first runs on `h1`: the iteration where `Cd`
becomes flat in upwind, and the iterations that van Leer needs.

### 5.4 Stop criterion

`system/FOs/FOrunTimeControl.j2`, the same file for both templates. For each
field in `coeffs_fields`, the `average` condition calculates a moving
average over `coeffs_window` iterations. `timeStart minIter` starts it at
`minIter`; in the central template, `controlDict_final` also includes it
only from `minIter`. The run stops when, for all the fields, the window is
complete and |average − current value| < `coeffs_variation`. `endIter` is
the second limit.

To change the tolerance of a running case, edit `tolerance` in its
`system/FOs/FOrunTimeControl` and `touch system/controlDict`. The solver
reads the file again and starts the average from zero, so the run cannot
stop before one more `coeffs_window`.

## 6. Turbulence

File: `constant/turbulenceProperties`.

```
simulationType  RAS;
RASModel        kOmegaSST;
turbulence      on;
printCoeffs     on;
```

**Closure.** The Boussinesq hypothesis gives the Reynolds stress from an
eddy viscosity and the strain rate:

```
−⟨u'ᵢu'ⱼ⟩ = νt (∂Uᵢ/∂xⱼ + ∂Uⱼ/∂xᵢ) − (2/3) k δᵢⱼ
νt = a₁ k / max(a₁ ω, S F₂)
```

The denominator is the shear-stress limiter. It limits the production of νt
in adverse pressure gradients. These gradients occur on the boattail and at
the fin-body junction.

**References.**

- Menter (1994), "Two-Equation Eddy-Viscosity Turbulence Models for
  Engineering Applications", *AIAA Journal* **32**(8), 1598–1605.
- Menter, Kuntz and Langtry (2003), "Ten Years of Industrial Experience with
  the SST Turbulence Model", *Turbulence, Heat and Mass Transfer* **4**,
  625–632. OpenFOAM uses this version. When you compare with the literature,
  find which version the source uses.
- Source code: `src/TurbulenceModels/turbulenceModels/RAS/kOmegaSST/`.
  Default constants: `alphaK1 0.85`, `alphaK2 1.0`, `alphaOmega1 0.5`,
  `alphaOmega2 0.856`, `beta1 0.075`, `beta2 0.0828`, `betaStar 0.09`,
  `gamma1 5/9`, `gamma2 0.44`, `a1 0.31`, `b1 1.0`, `c1 10`, `F3 no`.

**Wall treatment** on `cone`, `walls`, `tail`, `fins`:

| Field | Condition | Reason |
|---|---|---|
| `nut` | `nutUSpaldingWallFunction` | Spalding's law is continuous in y+. The first cell can be in the viscous sublayer, the buffer layer or the log layer. Spalding (1961), *J. Appl. Mech.* **28**(3), 455–458. |
| `k` | `kLowReWallFunction` | Also continuous in y+. |
| `omega` | `omegaWallFunction` | Blends the sublayer and log-layer values of ω. |
| `alphat` | `compressible::alphatWallFunction`, `Prt 0.85` | Turbulent thermal diffusivity from νt. |

These conditions need wall patches of type `wall`. Refer to MESH.md §6.1.

The mesh does not set y+. The velocity changes by a factor of 6 from M 0.3
to M 1.8. Measure y+ on each case. `Allrun` shows the command that reads
y+ from the log.

## 7. Thermophysical properties

File: `constant/thermophysicalProperties`.

| Entry | Value | Function |
|---|---|---|
| `type` | `hePsiThermo` | Energy-based thermodynamics with ψ = 1/(RT). `rhoCentralFoam` accepts only `psiThermo`. |
| `mixture` | `pureMixture` | One species. |
| `transport` | `sutherland` | μ(T) = As √T / (1 + Ts/T). |
| `thermo` | `hConst` | Constant Cp. h = Cp (T − Tstd) + Hf. |
| `equationOfState` | `perfectGas` | p = ρ R T. |
| `energy` | `sensibleInternalEnergy` | The solver solves e. `rhoCentralFoam` requires it. |

| Constant | Value | Calculation |
|---|---|---|
| `molWeight` | 28.96 kg/kmol | R = 8314.46 / 28.96 = 287.05 J/(kg K) |
| `Cp` | 1005 J/(kg K) | γ = Cp / (Cp − R) = 1005 / 717.95 = 1.400 |
| `As` | 1.4792e-6 kg/(m s √K) | At 288.15 K: μ = 1.790e-5 Pa s |
| `Ts` | 116 K | |

White, *Viscous Fluid Flow*, 3rd ed., eq. 1-36, gives Sutherland's law as
μ = μref (T/Tref)^1.5 (Tref + S)/(T + S). The OpenFOAM form
As √T / (1 + Ts/T) = As T^1.5 / (T + Ts) is the same equation.

**Limit of hConst.** The stagnation temperature is
T₀ = T (1 + (γ − 1)/2 · M²).

| M | T₀ | Cp error |
|---|---|---|
| 0.3 | 293 K | negligible |
| 0.8 | 325 K | < 0.5 % |
| 1.2 | 371 K | ~0.8 % |
| 1.8 | 475 K | ~1.8 % |
| 2.5 | 648 K | ~4 % |

Above M 2.5, replace `hConst` with `janaf`. The coefficients for air are in
`$FOAM_TUTORIALS/compressible/rhoCentralFoam/biconic25-55Run35/constant/thermophysicalProperties`.

`Prt = 0.85` is the turbulent Prandtl number for air. Kays (1994),
*J. Heat Transfer* **116**, 284–295.

## 8. Temperature limits

There are none in `case-central/`. `rhoCentralFoam` in v2412 does not read
`fvOptions`, so a `limitTemperature` entry has no effect. `case-pimple/` has
one (§16.4).

The initial transient can move T out of range. With `perfectGas`, a negative
T gives an imaginary speed of sound, and the solver stops with a floating
point exception. The Courant ramp (§5) is the only protection.

## 9. Boundary conditions

The number of characteristics that go in and go out on each boundary sets
the conditions. `U` and `T` use the same condition at all Mach numbers. `p`
does not. `0.orig/p` selects the branch with `#if ${{ $Minf > 1 }}`.

| Field | `inlet`, `outlet`, `box` | Function |
|---|---|---|
| `U` | `freestreamVelocity` | Fixed where the flow goes in. Extrapolated where it goes out. |
| `T`, `k`, `omega` | `inletOutlet` | Same. |
| `nut`, `alphat` | `calculated` | From the turbulence model. |

**p, M ≤ 1.** `freestreamPressure` on `inlet`, `outlet` and `box`. It
extrapolates where the flow goes in and is fixed where the flow goes out. In
subsonic flow, one acoustic wave goes in through the outlet.

**p, M > 1.**

| Boundary | Condition | Reason |
|---|---|---|
| `inlet` | `fixedValue` | All characteristics go in. |
| `outlet` | `zeroGradient` | All characteristics go out. |
| `box` | `waveTransmissive` | The normal velocity is approximately 0. Waves cross the boundary in the two directions. The condition moves p out at u + c and relaxes it to `fieldInf` over `lInf = 10` m. |

In v2412, `freestreamPressure` has no supersonic mode. On a supersonic
outlet it fixes p. This is the reason for the branch.

`waveTransmissive` is the non-reflecting condition of Poinsot and Lele
(1992), *J. Comput. Phys.* **101**, 104–129. At M 1.8 the Mach cone crosses
`box` (MESH.md §5). A `fixedValue` there reflects the wave into the domain.

**Walls.**

| Field | Condition |
|---|---|
| `U` | `noSlip` |
| `p` | `zeroGradient` |
| `T` | `zeroGradient` (adiabatic) |

The wall is adiabatic. At M 1.8 the recovery temperature is approximately
470 K. In a 4 s motor burn, the airframe temperature does not get to this
value. Zero heat flux is thus the correct limit.

The `symm` entries in `0.orig/` apply only to `half` and `quarter` meshes.
On a `full` mesh they have no effect.

## 10. Function objects

File: `system/FOs/FOcommon`, included by both templates. `residuals` is in
`controlDictBase.j2` and, for `rhoPimpleFoam`, in `case-pimple/system/controlDict.j2`,
because the solved fields differ.

| Name | Type | Output | Interval |
|---|---|---|---|
| `forceCoeffs1` | `forceCoeffs` | `Cd`, `Cl`, `Cm` and others | 50 time steps |
| `forces1` | `forces` | Forces and moments, pressure and viscous parts | 50 time steps |
| `MachNo` | `MachNo` | Mach number field | write time |
| `yPlus` | `yPlus` | y+ on the wall patches | write time |
| `residuals` | `solverInfo` | Residuals of `U`, `e`, `k`, `omega` (and `p` in pimple) | 50 time steps |
| `wallShearStress` | `wallShearStress` | Wall shear stress on `wallPatches` | write time |
| `Cp` | `pressure`, `staticCoeff` | Static pressure coefficient | write time |
| `surfaces` | `surfaces`, VTK | p, U, T, ρ, Mach, Cp, τw, y+ on `symm` and on the wall patches | end of the run |

**Reference values of the force coefficients.**

| Entry | Value | Source |
|---|---|---|
| `patches` | `$wallPatches` | `constant/meshInfo` |
| `rho` | `rho` | The density field. p is in Pa. |
| `rhoInf` | `$rhoInf` | `flowDerived`. Used for the dynamic pressure. |
| `magUInf` | `$Uinf` | `flowDerived` |
| `Aref` | `$Aref` | `constant/meshInfo` |
| `lRef` | `$lRef` | `constant/meshInfo` |
| `CofR` | `(0 0 0)` | The nose tip. MESH.md §3. |
| `pitchAxis` | `(0 0 1)` | z axis |
| `dragDir`, `liftDir` | from α and β | `flowDerived` |

**Time control.**

| Entry | Value | Function |
|---|---|---|
| `startFrom` | `latestTime` | A new `./Allrun` continues the run. |
| `endTime` | `endIter` | Iteration limit. §5.3. |
| `deltaT` | 1 | With `localEuler` the time counts iterations. §5.1. |
| `adjustTimeStep` | no | Not used with `localEuler`. |
| `maxDeltaT` | 1 s | Maximum local time step: no limit. |
| `writeInterval` | 500 time steps | Interval of the written times. |
| `purgeWrite` | 1 | Keeps the last written time. |
| `writeFormat` | binary | |

## 11. Discretization schemes

`system/fvSchemesVanLeer` contains the full set. `system/fvSchemes` and
`system/fvSchemesUpwind` are the same file: they include `fvSchemesVanLeer`
and set the three `reconstruct` entries to `upwind`. The ramp copies
`fvSchemesUpwind` at the start and `fvSchemesVanLeer` at iteration
`vanLeer`. Do not edit `system/fvSchemes`.

| Term | Scheme | Reason |
|---|---|---|
| `fluxScheme` | `Kurganov` | §4 |
| `reconstruct(rho)`, `reconstruct(T)` | `vanLeer` | §4 |
| `reconstruct(U)` | `vanLeerV` | §4 |
| `ddt` | `localEuler` | Local time step, pseudo-steady. §5.1. |
| `grad` | `Gauss linear` | |
| `div(tauMC)` | `Gauss linear` | Deviatoric part of the viscous stress. |
| `div(phi,k)`, `div(phi,omega)` | `Gauss upwind` | k and ω must stay positive. |
| `laplacian` | `Gauss linear limited corrected 0.33` | Non-orthogonal correction, limited. |
| `snGrad` | `limited corrected 0.33` | Same. |
| `wallDist` | `meshWave` | Wall distance for the SST model. |

## 12. Linear solvers

File: `system/fvSolution`.

| Fields | Solver | Tolerance | `relTol` |
|---|---|---|---|
| `rho`, `rhoU`, `rhoE` | `diagonal` | – | – |
| `U`, `e` | `smoothSolver`, `symGaussSeidel`, 2 sweeps | 1e-9 | 0.01 |
| `k`, `omega` | `smoothSolver`, `symGaussSeidel`, 1 sweep | 1e-9 | 0.1 |

The conserved variables are explicit. `diagonal` divides by the diagonal.
Only the implicit diffusive corrections and the turbulence equations use a
real linear solver.

## 13. Parallel decomposition

File: `system/decomposeParDict`. `numberOfSubdomains 8`, `method scotch`.
`newCase.sh --np` and `./Allrun <N>` change `numberOfSubdomains`.

## 14. Origin of the configuration

The template follows
`$FOAM_TUTORIALS/compressible/rhoCentralFoam/biconic25-55Run35` and the
Courant ramp of the OpenFOAM ToolChain. Differences from the tutorial:

| Item | Tutorial | Template | Reason |
|---|---|---|---|
| Turbulence | laminar | `kOmegaSST` | The rocket flow is turbulent. |
| Wall conditions | `maxwellSlipU`, `smoluchowskiJumpT` | `noSlip`, adiabatic | The tutorial is a low-density flow. The rocket flies at sea-level density. |
| `laplacian`, `snGrad` | not limited | `limited corrected 0.33` | Non-orthogonal cells. |
| `ddt` | `localEuler` | `localEuler` | Until 2026-10-12 the template used `Euler`, because a test on the previous mesh printed Courant numbers of 1e6 and diverged. The Courant number is not meaningful with `localEuler` (§5.1); that test had no ramp and had negative cell volumes. On `h1` at M 1.8 with the ramp, `localEuler` is stable (§16.6). |

## 15. Open items

| Item | Status |
|---|---|
| Validation | Not done. VALIDATION.md. |
| Wall patch types | `cartesianMesh` sets `wall`; `Allmesh` checks it. MESH.md §6.1. |
| Stop criterion | Window of 1000 iterations, from the ToolChain. Adjust it from the `h1` runs. §5.3. |
| Schedule of the ramp | `vanLeer`, `minIter` from the ToolChain. Adjust them from the `h1` runs. §5.3. |
| Cold cell at the base | On `h1` at M 0.3, one cell near x = 2.925 m, r = 0.06 m goes to the 100 K limit of `limitTemperature`. §16.6. |
| Motor jet | The base is a wall. A jet needs a `nozzle` patch and a total-pressure and total-temperature inlet. The base pressure with the motor on is different. |
| y+ | Not set by the mesh. Measure it on each case. |

## 16. rhoPimpleFoam template (M < 1)

### 16.1 Files

`case-pimple/` holds only the files that differ from `case-central/`.
`newCase.sh --solver pimple` copies `case-central/` without the files in
`case-pimple/central.exclude`, then copies `case-pimple/` on top. The other
files, `0.orig/`, `constant/`, `flowConditions`, `flowDerived`,
`decomposeParDict`, `FOs/FOcommon` and `FOs/FOrunTimeControl.j2`, are the
ones of `case-central/`: change them there.

| File | Function |
|---|---|
| `config.json` | `maxCo`, `minIter`, `endIter`, stop criterion. §5.3. |
| `system/controlDict.j2` | Time in iterations, function objects, `residuals` with `p`. |
| `system/fvSchemes` | §16.2. |
| `system/fvSolution.j2` | §16.3. |
| `system/fvOptions` | §16.4. |

### 16.2 Schemes

From `templates/rhoPimpleFoamTransonic` of the ToolChain, with the
`laplacian` and `snGrad` of `case-central/` for the non-orthogonal cells.

| Term | Scheme | Reason |
|---|---|---|
| `ddt` | `localEuler` | §5.1 |
| `grad` | `cellLimited Gauss linear 1`; 0.333 for `p` and `e` | Bounded gradients for the second-order convection. |
| `div(phi,U)` | `Gauss vanLeerV` | TVD, second order away from the shocks of M 0.8. |
| `div(phi,e)`, `div(phid,p)` | `Gauss limitedLinear 1` | Bounded. `phid` is the convective flux of p in the `transonic` formulation. |
| `div(phi,K)` | `Gauss linear` | Kinetic energy. |
| `div(phiv,p)` | `Gauss upwind` | Pressure work in the e equation. |
| `div(phi,k)`, `div(phi,omega)` | `Gauss upwind` | k and ω must stay positive. Same as central. |
| `laplacian`, `snGrad` | `limited corrected 0.33` | Same as central. |

### 16.3 PIMPLE

| Entry | Value | Reason |
|---|---|---|
| `transonic` | yes | Density-weighted pressure equation, needed for the compressibility at M 0.6–0.8. Valid at M 0.3. |
| `nOuterCorrectors` | 1 | Pseudo-steady march: one PISO pass per iteration. |
| `nCorrectors` | 2 | |
| `nNonOrthogonalCorrectors` | 1 | Maximum non-orthogonality 80°. |
| `pMinFactor`, `pMaxFactor` | 0.05, 6 | Clip limits of p, relative to the boundary values. 6 covers the stagnation pressure behind a normal shock at M 1.8. |
| `maxCo`, `rDeltaTSmoothingCoeff`, `rDeltaTDampingCoeff`, `maxDeltaT` | 0.5, 0.1, 0.9, 1 | §5.1. The ToolChain writes these in `controlDict`, where `rhoPimpleFoam` does not read them. |

Linear solvers: `p` with `GAMG`, tolerance 1e-8, `relTol` 0.01; `rho`, `U`,
`e`, `k`, `omega` with `smoothSolver` `symGaussSeidel`, `relTol` 0.1.

### 16.4 Temperature limit

`limitTemperature` between 100 K and 1000 K. `rhoPimpleFoam` reads
`fvOptions`. The free stream is at 236–276 K and the stagnation
temperature at M 1.8 is 475 K. The limit only protects the start: count the
cells that it limits in the log (`LimitedCells`).

### 16.5 Turbulence, gas, boundary conditions

The same files as `case-central/`: §6, §7, §9. The `p` branch at
`Minf > 1` also works with `rhoPimpleFoam`, but the sweep uses it only at
M < 1.

### 16.6 Tests on h1 (2026-10-12)

`h1`: half mesh, 2 198 939 cells, 8 processes on a Ryzen 7 250.

| Case | Result |
|---|---|
| Pimple, M 0.3, 150 iterations from the uniform free stream | Stable. 4.6 s per iteration, 90 s of start. `Cd` 0.61 at 100, 0.42 at 150, not converged. y+ mean (max): cone 6.6 (8.8), walls 3.1 (7.8), tail 2.1 (6.1), fins 11.1 (30.1). |
| Pimple, M 0.3, `--init` from the case above, 50 iterations | `mapFields` and the run work. `mapFields` warns about tet decomposition on about 100 faces of the layer cells. It also maps `rDeltaT`. |
| Pimple, M 0.3, 100 more iterations: original schemes, `div(phi,K)` upwind, `div(phi,e)` and `div(phi,K)` upwind | The same within 1 %: minimum of T 129 → 175 K, maximum 341 → 271 K, `Cd` 0.300. The schemes do not change the cold cell. |
| Pimple, M 1.8, 800 iterations from the uniform free stream | §17. `Cd` 0.69 (50), 0.54 (300), 0.47 (500 to 800). |
| Central, M 1.8, from the uniform free stream, three variants | Crash at iteration 74 to 256. §17. |
| Central, M 1.8, `--init` from the pimple run at 800, 400 iterations through the whole ramp | Stable. Minimum of T 136–218 K at the fin tip and fin root, maximum 481 K. `Cd` 0.62 in upwind, 0.51–0.53 in van Leer, not converged. |

At M 0.3, at iteration 1, `limitTemperature` limits 160 cells. After that,
0 to 2 cells. The unlimited minimum of T is in a layer cell near
x = 2.925 m, r = 0.06 m, between the fin trailing edge and the boattail
(aspect ratio 28, determinant 2e-4), almost at rest (U ≈ 5 m/s) with p near
the stagnation value. It rises slowly with the iterations, with any of the
schemes above: it is the start transient in a region with almost no flow.

### 16.7 Two solvers at the same Mach number

The sweep runs M 1.2 and M 1.8 on `h0` with both solvers (§17), so it gives
the difference between the solvers at those Mach numbers, on the coarse
mesh. On `h1` at M 1.8, after the tests above: pimple 0.47, central 0.52,
neither converged. To measure the step on a mesh of the study, add a pimple
case on `h1` started from `m18-h0-pimple`, or a central case at M 0.8.

## 17. Start of the supersonic cases

From the uniform free stream, the flow behind the base starts at
`Uinf`: the base acts as a piston and empties the region behind it. The
cells there are thin layer cells with a small local time step, and they
fill slowly. On `h1` at M 1.8, with the pimple template:

| Iteration | p at the base [Pa] | T at the base [K] |
|---|---|---|
| 305 | 7 900 | 166 |
| 455 | 57 000 | ~320 |
| 755 | 46 000 | 350–414 |

p∞ is 80 127 Pa. At the end p_base/p∞ ≈ 0.57, in the range of a power-off
base at M 1.8, and T at the base goes to the recovery temperature of the
adiabatic wall. During the first 300 iterations, T at the base stays at the
100 K limit of `limitTemperature`.

`rhoCentralFoam` reads no `fvOptions`, so it has no such limit. In three
tests from the uniform free stream (default `rDeltaTSmoothingCoeff` 0.1;
0.02; `maxDeltaT 1e-6`), T at the base fell to 0 and the solver stopped with
a floating point exception in the thermophysical model, at iterations 74 to
256. The parameters of the local time step do not change this.

Thus, for M > 1, the sweep runs the pimple template first on `h0` and
starts the central case from it with `--init` (WORKFLOW.md §5.6). `h1`, `h2`
and `h3` start from the central solution of the coarser mesh. Tested on
`h1`: pimple → central on the same mesh. Not tested yet: pimple on `h0`,
and central from a central solution on another mesh.
