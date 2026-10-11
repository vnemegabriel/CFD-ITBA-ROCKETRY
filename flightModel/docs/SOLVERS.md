# Solver and models

This document gives the dictionaries of the case template `case-central/`,
the models that they select, and the references. Use it to calculate the
case values by hand and to compare them with the case.

> **CAUTION** The template is not validated. Do not use the results before
> you do the procedures in [VALIDATION.md](VALIDATION.md).

## 1. Summary

| Item | Value |
|---|---|
| Template | `case-central/` |
| Solver | `rhoCentralFoam` |
| Formulation | Compressible, density-based, explicit |
| Range | M 0.3 to M 1.8, one template |
| Time | `Euler`, global time step from `maxCo`, Courant ramp (§5) |
| Fields in `0.orig/` | `U`, `p`, `T`, `k`, `omega`, `nut`, `alphat` |
| Turbulence | k-ω SST with wall functions |
| Gas | Air, perfect gas, constant Cp, Sutherland viscosity |
| Input | `Minf`, `pInf`, `Tinf`, `alpha`, `beta`, `Ti`, `nuRatio` |
| Initial condition | Uniform free stream |

One solver for all Mach numbers. With three solvers, the discretization
changes between flow regimes. A step in the Cd–Mach curve is then not
separable from the physics.

## 2. Files of the template

| File | Edit | Function |
|---|---|---|
| `system/flowConditions` | yes | Flight condition. §3. |
| `config.json` | yes | Parameters of the Courant ramp and of the stop criterion. §5. |
| `system/flowDerived` | no | Values calculated from `flowConditions`. §3. |
| `0.orig/*` | no | Initial and boundary conditions. §9. |
| `constant/thermophysicalProperties` | no | Gas model. §7. |
| `constant/turbulenceProperties` | no | Turbulence model. §6. |
| `constant/meshInfo` | no | `newCase.sh` replaces it with the `meshInfo` of the mesh. MESH.md §7. |
| `system/controlDictBase.j2` | yes, with care | Time control and function objects. §5, §10. |
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

## 5. Courant ramp and stop criterion

### 5.1 Function

A uniform free stream on a no-slip wall causes a large transient. The ramp
starts with a small `maxCo` (0.01) and first-order schemes. Then it increases
`maxCo` step by step. The method comes from the OpenFOAM ToolChain
(`templates/rhoCentralFoam`).

The `timeActivatedFileUpdate` function objects in `system/fileUpdater` do the
ramp. At each time step they compare the time with a schedule. At each
switch time they copy a file on `system/controlDict` or on
`system/fvSchemes`. `runTimeModifiable yes` makes the solver read the new
file.

### 5.2 config.json

`newCase.sh` renders all `*.j2` files with these values.

| Key | Default | Function |
|---|---|---|
| `tauCo` | `3.3e-8` | Time step in s that gives Co = 1 in the worst cell. |
| `vanLeer` | `2000` | Iteration of the change from upwind to van Leer. `-1` keeps upwind. |
| `minIter` | `3000` | Iteration of the change to `controlDict_final`. |
| `final_Co` | `0.5` | Maximum `maxCo` of all the ramp. |
| `coeffs_fields` | `["Cd"]` | Coefficients in the stop criterion. |
| `coeffs_variation` | `[0.0001]` | Absolute tolerance of each coefficient. |
| `coeffs_window` | `[2]` | Window of each coefficient, in flow-through times. |

At α ≠ 0 you can add `Cl` and `CmPitch` to `coeffs_fields`. The tolerance
is absolute, so it also applies to a coefficient with an average near 0.

### 5.3 Schedule

The schedule is in iterations:

| Phase | Iterations | `maxCo` |
|---|---|---|
| Upwind | 0, 20, 35, 50, 65, 80, 100 | 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.5 |
| van Leer | `vanLeer` + 0, 20, 35, 150, 300, 800, 2000 | 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.5 |
| Final | `minIter` | `final_Co`, with `runTimeControl` |

Each `maxCo` is limited to `final_Co`. The template removes the events after
`minIter`.

The case uses `Euler` in physical time, not `localEuler`. The template thus
changes each iteration into a time:

```
t = Σ nᵢ · maxCoᵢ · tauCo
```

The initial time step is `deltaT = 0.01 · tauCo / 1.2`. `rhoCentralFoam`
makes its first step at 1.2·`deltaT` and does not examine the Courant
number.

### 5.4 Measure tauCo

`tauCo` depends on the mesh and a little on the Mach number. Values on the
base mesh of `Allmesh` (2 144 494 cells), measured on the first steps from
the uniform free stream:

| M | `tauCo` |
|---|---|
| 0.3 | 3.31e-8 s |
| 1.8 | 2.71e-8 s |

`config.json` has `3.3e-8`. A `tauCo` larger than the real value makes each
phase of the ramp last more iterations than the schedule, which is the safe
direction. `maxCo` still limits the Courant number.

Measure `tauCo` again when you change the mesh:

1. Make a case and run some time steps.
2. In `log.rhoCentralFoam`, read `deltaT` and the max Courant number on the
   next line. The solver calculates the Courant number with that `deltaT`.
3. Calculate `tauCo = deltaT / Co_max`. All the steps from the uniform
   free stream give the same value.
4. Write the largest value in `case-central/config.json`.
5. Make the cases again with `newCase.sh`.

### 5.5 Stop criterion

`controlDict_final` adds `runTimeControl1`. For each field in
`coeffs_fields`, the `average` condition calculates a moving average with a
window of `coeffs_window` flow-through times:

```
window = coeffs_window · lBody / Uinf
```

The window is in physical time. The condition starts when the run changes to
`controlDict_final`, at `minIter`. The run stops when both are true:

1. The physical time since `minIter` is more than the window.
2. |average − current value| < `coeffs_variation`, for all the fields.

`endTime = 10 · lBody / Uinf` is the second limit:

| M | Flow-through time | `endTime` |
|---|---|---|
| 0.3 | 28.9 ms | 0.289 s |
| 0.8 | 10.9 ms | 0.109 s |
| 1.8 | 4.82 ms | 0.048 s |

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

There are none. `rhoCentralFoam` in v2412 does not read `fvOptions`, so a
`limitTemperature` entry has no effect.

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

File: `system/controlDictBase.j2`.

| Name | Type | Output | Interval |
|---|---|---|---|
| `forceCoeffs1` | `forceCoeffs` | `Cd`, `Cl`, `Cm` and others | 50 time steps |
| `forces1` | `forces` | Forces and moments, pressure and viscous parts | 50 time steps |
| `MachNo` | `MachNo` | Mach number field | write time |
| `yPlus` | `yPlus` | y+ on the wall patches | write time |
| `residuals` | `solverInfo` | Residuals of `U`, `e`, `k`, `omega` | each time step |
| `wallShearStress` | `wallShearStress` | Wall shear stress on `wallPatches` | write time |
| `Cp` | `pressure`, `staticCoeff` | Static pressure coefficient | write time |

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
| `endTime` | 10 · `lBody` / `Uinf` | Physical time limit: 10 flow-through times. |
| `deltaT` | `0.01·tauCo/1.2` | Initial time step. §5.3. |
| `adjustTimeStep` | yes | The time step follows `maxCo`. |
| `maxDeltaT` | 1e-5 s | Maximum time step. |
| `writeInterval` | 0.005 s, adjustable | Interval of the written times. |
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
| `ddt` | `Euler` | Global time step. `localEuler` did not limit the step on the previous mesh. §14. |
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
| `ddt` | `localEuler` | `Euler` | On the previous mesh, `localEuler` gave Courant numbers of 1e6 and diverged in less than ten iterations. That test had no ramp and had negative cell volumes. Do the test again on the cfMesh mesh. |

## 15. Open items

| Item | Status |
|---|---|
| Validation | Not done. VALIDATION.md. |
| Wall patch types | `cartesianMesh` sets `wall`; `Allmesh` checks it. MESH.md §6.1. |
| `tauCo` | Measured on the base mesh, from the uniform free stream. §5.4. |
| Stop criterion | Window of 2 flow-through times. §5.5. |
| Local time step | Not tested on the cfMesh mesh. §14. |
| Motor jet | The base is a wall. A jet needs a `nozzle` patch and a total-pressure and total-temperature inlet. The base pressure with the motor on is different. |
| y+ | Not set by the mesh. Measure it on each case. |
