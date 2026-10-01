# Supersonic Fin Flutter Solver

High-fidelity aeroelastic analysis of a composite CFRP rocket fin using Q4 Mindlin shell FEM, 2nd-order piston theory, and the p-k flutter method. Developed for IREC 2026 (Team 207, ITBA).

---

## What It Does

Given a composite laminate design and a rocket flight envelope, this solver:

1. Builds a finite element model of the fin (Q4 Mindlin shell, anisotropic CLT bending)
2. Extracts natural frequencies and mode shapes
3. Computes generalized aerodynamic forces (GAF) via 2nd-order piston theory
4. Solves the p-k eigenvalue problem across all supersonic flight conditions
5. Identifies the flutter speed, divergence speed, and safety margins

**Key results for the AR1 / β=5° / 8 mm laminate:**

| Quantity | Value | Margin vs flight |
|---|---|---|
| Natural frequencies | 237, 592, 1342, 1439, 2057, 2247 Hz | — |
| Divergence speed | 2850 m/s | 4.6× |
| Flutter speed (g=0.02) | 568 m/s | 0.92× |
| Flutter mechanism | Single-mode anti-damping (NOT bending-torsion) | — |

Bending-torsion flutter is **fully suppressed** by the D16 aeroelastic tailoring (β=5° fiber rotation). The remaining single-mode flutter speed scales as V_f ≈ 4334·√g m/s and depends critically on the structural loss factor g.

---

## Repository Layout

```
fin_flutter/
├── README.md
├── LICENSE
├── IREC2026_improvement_plan.md
├── Recursos/                          # Reference PDFs (Dowell, Bisplinghoff, Jones, etc.)
└── v2/
    └── supersonic_fin_flutter_matlab/
        ├── mainFlutterSolver.m        ← top-level pipeline; run from this directory
        ├── CLAUDE.md                  ← detailed technical context and design notes
        ├── data/
        │   ├── lam.json               ← CLT output: D-matrix, density, thickness (β=5°)
        │   └── flight_data.csv        ← telemetry: time, altitude_ft, Vz_ms (229 points)
        ├── functions/
        │   ├── isaAtmosphere.m        ← ISA atmosphere model
        │   ├── GenerarMallaAleta.m    ← Q4 swept trapezoidal mesh generator
        │   ├── CalcularRigidezQLLL.m  ← Q4 Mindlin element stiffness (anisotropic CLT)
        │   ├── assembleGlobalStiffness.m
        │   ├── assembleGlobalMass.m
        │   ├── applyDirichletBCs.m    ← full root clamp
        │   ├── modalAnalysis.m        ← eigensolver, mass-normalized modes
        │   ├── pistonTheoryGAF.m      ← 2nd-order piston theory GAF matrix
        │   └── pkSolveFlutter.m       ← p-k flutter + divergence solver
        ├── python_core/
        │   └── inplaneG_v5.py         ← CLT laminate generator → produces lam.json
        └── results/                   ← auto-generated plots and flutter.mat
```

---

## Quick Start

### Requirements

- MATLAB R2023b or later (Control System Toolbox for eigenvalue solvers)
- Python 3.8+ (only needed to regenerate `lam.json`)

### Run the Flutter Solver

```matlab
% In MATLAB, from the v2/supersonic_fin_flutter_matlab/ directory:
mainFlutterSolver
```

The script loads `data/lam.json` and `data/flight_data.csv`, runs the full pipeline, and writes plots and `results/flutter.mat`.

### Regenerate the Laminate (optional)

```bash
cd v2/supersonic_fin_flutter_matlab/python_core
python inplaneG_v5.py --beta 5 --thickness 8.0 --json ../data/lam.json
```

---

## Analysis Pipeline

### Step 1 — Laminate Design (CLT)

`inplaneG_v5.py` builds the ABD stiffness matrix for the AR1 layup (DB300 ±45° NCF + GA90R 0/90° woven, T700/Epoxy). A global fiber rotation by β degrees introduces a non-zero D16 bending-torsion coupling term.

Current laminate (`lam.json`, β=5°, t=8 mm):

| Property | Value | Units |
|---|---|---|
| D11 | 2018 | N·m |
| D22 | 2040 | N·m |
| D66 | 657  | N·m |
| D16 | +41.3 | N·m |
| Density | 1501 | kg/m³ |

### Step 2 — FEM Model

`GenerarMallaAleta` builds a Q4 swept trapezoidal mesh (default 24×12 elements) matching the fin geometry. `CalcularRigidezQLLL` assembles the element stiffness matrix using the anisotropic CLT D-matrix directly as the constitutive law (Mindlin shell, 6 DOF/node, selective reduced integration).

**Fin geometry:**

| Parameter | Value |
|---|---|
| Root chord (cr) | 0.300 m |
| Tip chord (ct)  | 0.150 m |
| Semi-span        | 0.160 m |
| LE sweep angle   | 57.4°   |

Root nodes are fully clamped (all 6 DOF fixed).

### Step 3 — Modal Analysis

`modalAnalysis` solves the generalized eigenvalue problem K·Φ = M·Φ·Ω² and mass-normalizes the mode shapes (Φᵀ·M·Φ = I). Six modes are retained.

### Step 4 — Generalized Aerodynamic Forces

`pistonTheoryGAF` computes the complex GAF matrix Q(k) for a set of reduced frequencies using the 2nd-order Lighthill (1953) piston theory pressure kernel:

```
p_j = (2·q_∞/β) · (i·k/b_ref · w_j + ∂w_j/∂x)
Q_ij = ∫ φ_i · p_j · dA   [2×2 Gauss quadrature per element]
```

where β = √(M²−1), b_ref = MAC/2 = 0.1167 m, and k = ω·b/U is the reduced frequency. Q is intentionally **non-symmetric** — the anti-symmetric part encodes the bending-torsion coupling responsible for classical flutter.

### Step 5 — p-k Flutter Solution

`pkSolveFlutter` sweeps dynamic pressure q from zero to beyond divergence. At each q step it solves the complex eigenvalue problem:

```
A(p,q) = q·(Q0 + i·k·Q1) − Ω²_c
p = sqrt(eig(A)),   Im(p) ≥ 0
```

where Ω²_c = diag(ωₙ²·(1 + i·g)) is the complex modal stiffness (g = structural loss factor). Flutter occurs when Re(p) first crosses zero. Divergence is extracted separately from the quasi-steady (k=0) eigenvalue problem.

**Key solver parameter:**

| Parameter | Location | Default | Effect |
|---|---|---|---|
| `g_struct` | `pkSolveFlutter.m` | 0.02 | V_flutter ∝ √g — dominant sensitivity |
| `nx, ny` | `mainFlutterSolver.m` | 24, 12 | Mesh resolution |
| `nModes` | `mainFlutterSolver.m` | 6 | Modal truncation |

---

## Aeroelastic Tailoring

The β=5° fiber rotation introduces D16 = +41.3 N·m. For an aft-swept fin (Λ=57.4°), this creates a **wash-out** mechanism: bending deflection induces a nose-down torsional rotation, reducing the effective angle of attack and opposing further deformation.

Physical confirmation from the GAF matrix: Q0(1,2)·Q0(2,1) = −0.24 < 0 (stabilizing cross-coupling). Classical bending-torsion flutter requires this product to be positive — it is not, confirming suppression.

The remaining flutter mode is **single-mode aerodynamic anti-damping**, driven by the ẇ/U pressure term in piston theory. Its flutter speed has a closed-form approximation:

```
V_flutter ≈ 4334·√g   m/s   (at worst-case Mach 1.838)
```

| g (structural loss factor) | V_flutter | Margin vs V_flight = 618 m/s |
|---|---|---|
| 0.005 | 281 m/s | 0.45× |
| 0.010 | 397 m/s | 0.64× |
| 0.020 | 562 m/s | 0.91× |
| 0.030 | 688 m/s | 1.11× |
| 0.050 | 888 m/s | 1.44× |

The value g = 0.02 is a conservative lower bound for T700/epoxy from DMA literature. Experimental characterization is required for flight certification.

---

## Flight Envelope

`data/flight_data.csv` contains 229 supersonic telemetry points (Ma 1.05–1.85, 0–5 km altitude) covering the ascent max-q phase. The solver evaluates flutter margins at all 229 conditions. The worst-case point occurs at Mach 1.838, where the flight velocity reaches 618 m/s.

---

## Theory References

| Source | Used for |
|---|---|
| [Lighthill (1953) *J. Aeronautical Sciences* 20(6):402–406](https://arc.aiaa.org/doi/10.2514/8.2657) | 2nd-order piston theory pressure kernel |
| [Bisplinghoff, Ashley & Halfman (1955) *Aeroelasticity*, Dover ed.](https://store.doverpublications.com/products/9780486691893) | p-k method formulation |
| [Weisshaar (1981) *J. Aircraft* 18(8):669–676](https://doi.org/10.2514/3.57542) | Aeroelastic tailoring via D16/global β rotation |
| [Jones (1946) *NACA Report 835*](https://ntrs.nasa.gov/citations/19930091913) | Low-AR fin supersonic aerodynamic theory |
| Jones (1999) *Mechanics of Composite Materials* (2nd ed.) | CLT laminate theory |
| [Halpin & Tsai (1969) *AFML-TR-67-423*](https://apps.dtic.mil/sti/tr/pdf/ADA306357.pdf) | Composite micromechanics (Halpin-Tsai) |
| [MIL-A-8870C (1993)](http://everyspec.com/MIL-SPECS/MIL-SPECS-MIL-A/MIL-A-8870C_6746/) | Flutter and divergence certification requirements |

---

## Open Items

- **Structural damping:** g not yet measured for this laminate. DMA test needed for certification. V_flutter is highly sensitive to g (scales as √g).
- **Mesh convergence:** 24×12 mesh not formally verified. Recommended: refine to 48×24 and compare natural frequencies.
- **Mode 5 (2057 Hz):** Membrane-dominated mode with zero aerodynamic coupling — confirmed correct, not visualized explicitly.

---

## License

MIT License. See [LICENSE](LICENSE) for details.
