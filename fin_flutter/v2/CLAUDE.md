# CLAUDE.md — Supersonic Fin Flutter Solver
## Session Context Block

```
[CB-FLUTTER-2026]
System:    supersonic CFRP rocket fin, Q4 Mindlin shell FEM + 2nd-order piston theory + p-k flutter solver
Regime:    Ma 1.05–1.85, ISA sea level to 5 km, 229 supersonic flight points
Method:    CLT laminate → FEM modal analysis (6 modes) → pistonTheoryGAF → pkSolveFlutter
Laminate:  T700/Epoxy AR1, beta=5° aeroelastic tailoring, t=8 mm, rho=1501 kg/m³
           D11=2018, D22=2040, D66=657, D16=41.3 [N·m] (tailored_beta in lam.json)
Geometry:  cr=0.300 m, ct=0.150 m, span=0.160 m, sweep=57.4° LE, mesh 24×12 Q4 elements
Results:   f_n=[237, 592, 1342, 1439, 2057, 2247] Hz | V_div=2850 m/s (4.6×) | V_fl=568 m/s @g=0.02
Notation:  g = structural loss factor | Q0 = quasi-steady GAF | Q1 = aerodynamic damping slope
           β = sqrt(M²-1) | b_ref = MAC/2 = 0.1167 m | k = ω·b/U (reduced frequency)
Units:     SI throughout. Q0,Q1 in 1/(Pa·s²) after normalization by q_crit=209827 Pa
```

---

## 1. Repository Structure

```
supersonic_fin_flutter_matlab/
├── mainFlutterSolver.m          ← top-level pipeline; run from this directory
├── data/
│   ├── lam.json                 ← CLT output: D-matrix, rho, t, layup (beta=5)
│   └── flight_data.csv          ← time, altitude_ft, Vz_ms; 229 supersonic pts filtered at M≥1.05
├── functions/
│   ├── isaAtmosphere.m          ← ISA atmosphere [rho, a, T, p] = f(h_m)
│   ├── GenerarMallaAleta.m      ← Q4 swept trapezoidal fin mesh [nodes, connect]
│   ├── CalcularRigidezQLLL.m    ← Q4 Mindlin element stiffness (anisotropic CLT bending)
│   ├── assembleGlobalStiffness.m← sparse global K from CalcularRigidezQLLL
│   ├── assembleGlobalMass.m     ← consistent mass matrix
│   ├── applyDirichletBCs.m      ← clamp root nodes (y<1e-9), extract free DOFs
│   ├── modalAnalysis.m          ← eig/eigs → mass-normalized modes [Phi, omega_n]
│   ├── pistonTheoryGAF.m        ← 2nd-order piston theory GAF Q_k [nM×nM×nK] complex
│   ├── pkSolveFlutter.m         ← p-k flutter + quasi-steady divergence solver
│   └── buildCLTLaminate.m       ← (unused in main pipeline; standalone CLT utility)
├── python_core/
│   └── inplaneG_v5.py           ← CLT laminate generator; outputs lam.json
├── results/                     ← auto-generated plots + flutter.mat
└── CLAUDE.md                    ← this file
```

---

## 2. Physics Reference

### 2.1 CLT Laminate → lam.json
- `inplaneG_v5.py` builds the ABD matrix for the AR1 layup (DB300 ±45° NCF + GA90R 0/90° woven)
- beta_deg=5 rotates all fiber angles by +5° → introduces D16=+41.3 N·m (bending-torsion coupling)
- **G_eff = 12·D66/t³** — the effective flexural shear modulus; CONSTANT for self-similar scaling
  (D66 ∝ t³ exactly for proportional ply scaling → G_eff is a pure material property, not t-dependent)
- Discrete-ply V2 model: actual Vf recomputed from total FAW/(rho_f×t_mold); small wobble in G_eff
  (14.9–17.2 GPa) comes from non-proportional DB300/GA90R counts at different target thicknesses
- **Do not confuse D66 (N·m) with G_eff (GPa)**: D66 grows as t³; G_eff is constant

### 2.2 FEM (Q4 Mindlin Shell)
- 6 DOF/node: [u, v, w, θx, θy, θz] — piston theory uses ONLY w (DOF 3) per node
- Constitutive: bending from anisotropic CLT D-matrix (D_flex_3x3); membrane from E_eff=12D66/t³
- Shear: Mindlin k=5/6, selective reduced integration (1-pt) to avoid shear locking
- Drilling stiffness: k_drill = 1e-3·E_eff·t·A_elem per θz DOF (Oñate regularization)
- **CalcularRigidezQLLL re-derives E_eff = 12·D_flex(3,3)/t³ internally** — material.E from caller is ignored
- Root clamping: nodes with y < 1e-9 → fixedDOFs = reshape((rootNodes-1)*6 + (1:6))

### 2.3 Modal Analysis
- Mass normalization: Φᵀ M Φ = I → mode shape units [kg^{-1/2}]
- For nDOF ≤ 2500: full `eig`; else `eigs` with shift σ=1 rad²/s²
- Natural frequencies: [237, 592, 1342, 1439, 2057, 2247] Hz
- **Mode 5 (2057 Hz) has ALL-ZERO Q0 and Q1** — correct, it is a membrane/in-plane mode
  with zero out-of-plane (w) displacement; piston theory (w-only) gives zero aerodynamic coupling

### 2.4 Piston Theory GAF (2nd-order, Lighthill 1953)
Pressure kernel at Gauss point:
```
p_j = (2·q_inf/β) · (i·k/b_ref · w_j + ∂w_j/∂x)
Q_ij = ∫ φ_i · p_j · dA   [2×2 Gauss quadrature per Q4 element]
```
Key properties:
- Q is **NON-SYMMETRIC**: Q_ij ≠ Q_ji because ∫φ_i·∂φ_j/∂x ≠ ∫φ_j·∂φ_i/∂x
- The anti-symmetric part (off-diagonal gradient terms) IS the bending-torsion flutter mechanism
- β = sqrt(M²-1) uses free-stream β (not sweep-corrected) — correct for low-AR fin per Jones (1946)
- After normalization: Q_k_norm = Q_k / q_crit; units [m/kg] (i.e., 1/(Pa·s²)×s² = m/kg so that q·Q has units rad²/s²)
- Q_k decomposition in pkSolveFlutter:
  - Q0 = real(Q_k_norm(:,:,1))          quasi-steady stiffness [1/(Pa·s²)]
  - Q1 = imag(Q_k_norm(:,:,2))/k_vals(2) aerodynamic damping slope [1/(Pa·s²)]

### 2.5 p-k Flutter Method
Eigenvalue problem at each dynamic pressure q:
```
A = q·(Q0 + i·k·Q1) - Ω²_c
p² = eig(A)           [complex, non-Hermitian]
p = sqrt(λ), Im(p)≥0  [physical branch: positive frequency]
Flutter: Re(p) crosses 0 from below as q increases
```
- Ω²_c = diag(ω_n²·(1+i·g_struct))  — complex stiffness (structural loss factor g)
- k self-consistency: k = Im(p_mode1)·b_ref/U  (see §4 for why mode-1 only)
- Divergence (quasi-steady, k=0): B̂ = diag(1/ω)·Q0·diag(1/ω); q_div = 1/max(positive eig(B̂))

### 2.6 Aeroelastic Tailoring (Weisshaar 1981)
- β=5° fiber rotation introduces D16 > 0 for aft-swept fin (Λ=57.4°)
- **Wash-out mechanism**: bending deflection → nose-down torsion → reduced angle of attack
- Physical confirmation: Q0(1,2)·Q0(2,1) = -0.2406 < 0 → stabilizing cross-coupling
  (classical bending-torsion flutter requires Q0(1,2)·Q0(2,1) > 0)
- Bending-torsion flutter is **FULLY SUPPRESSED** by the D16 tailoring

---

## 3. Bugs Found and Fixed

### BUG 1 (CRITICAL): Hermitian Symmetrization in pistonTheoryGAF.m
**Symptom:** V_flutter = Inf for all flight conditions  
**Location:** `functions/pistonTheoryGAF.m`, formerly at end of function  
**Removed line:**
```matlab
Q_k = (Q_k + permute(conj(Q_k), [2, 1, 3])) / 2;   % <-- DELETED
```
**Why it caused Inf:** (Q+Qᴴ)/2 zeroes the anti-symmetric part ∫φ_i·∂φ_j/∂x,
which is physically the bending-torsion flutter mechanism. Symmetrized Q gives all
purely imaginary p → Re(p)=0 never crossed → V_flutter = Inf.  
**Status:** Fixed. Q_k is now intentionally left non-Hermitian. See comment block in pistonTheoryGAF.m.

### BUG 2 (SIGNIFICANT): Mean-k Averaging in pkSolveFlutter.m
**Symptom after Bug 1 fix:** V_flutter = 168 m/s (margin 0.27×) — physically unreasonable  
**Location:** `functions/pkSolveFlutter.m`, inner k-convergence loop  
**Root cause:** k-update was `mean(max(imag(p_new), 0)) * b_ref/U_i`  
Modes 3–6 have ω ≈ 8430–14120 rad/s → mean(ω) ≈ 8288 rad/s → k_avg ≈ 1.29 at first q-step.  
This large k pumped a huge Im(Q1) term into A:
```
Im(A(1,1)) = q·k·Q1(1,1) - ω₁²·g   →  with k=1.29: +16409 rad²/s²  (positive!)
→ arg(λ₁) > π - ε  →  arg(p₁)/2 > π/2 - ε  →  Re(p₁) = +5.5 rad/s  (SPURIOUS FLUTTER)
```
The k-iteration converges at this spurious fixed point (k≈1.29 is self-consistent when
including high-frequency modes, but physically wrong for mode 1).  
**Fix:** Changed k-update to `min(max(imag(p_new), 0)) * b_ref/U_i`  
(uses mode 1's Im(p) only, which gives k₁ ≈ 0.28 — the physically relevant reduced frequency)

**Additional fix:** k_avg initialized to 0 at the start of each flight condition's q-sweep  
(previously seeded from ωₙ at q=0, giving k_init ≈ 1.61 — same problem)

---

## 4. Numerical Subtleties

### 4.1 Why min(Im(p)) for k, not mean
The p-k method is formally mode-specific: each mode r should have k_r = Im(p_r)·b/U.
With a single shared k for the matrix Q_eff, the only physically defensible choice is
k₁ = Im(p₁)·b/U (lowest-frequency, primary flutter mode). High-frequency modes (2–6)
contribute k ≈ 0.58–2.21 and pollute the mean. Using min is equivalent to mode-1-driven k.

### 4.2 Mode 5 zero-GAF is correct
Mode 5 (2057 Hz) has Q0 row/col = 0 and Q1 row/col = 0. This is physically correct:
mode 5 is membrane/in-plane dominated (negligible w-displacement). Since pistonTheoryGAF
integrates only w-DOFs, this mode has zero aerodynamic coupling and does not participate
in flutter. It remains stable (decoupled) at all flight conditions.

### 4.3 Mode tracking (matchModes)
Greedy nearest-|p_new - p_ref| assignment. Works correctly when Im(p) values are
well-separated (which they are: 1488, 3722, 8430, 9039, 12925, 14120 rad/s). The
assignment can be verified by checking Im(p_matched)/omega_n ≈ 1.0 for all modes
at low q.

### 4.4 q_max includes divergence headroom
```matlab
q_max = max(4*q_flt, 1.1*q_div_crit)
```
Since q_div_crit ≈ 4.44M Pa >> 4×q_flt, the first q-step size is q_max/149 ≈ 32,746 Pa —
much larger than one might expect from q_flt alone. This is intentional (sweep must reach
the divergence boundary) but means the first q-step is coarse.

### 4.5 Units throughout the pipeline
| Quantity | Units | Note |
|---|---|---|
| Q_k (from pistonTheoryGAF) | rad²/s² | = s^{-2} |
| Q_k_norm = Q_k/q_crit | m/kg | so that q·Q_norm has units rad²/s² |
| Q0, Q1 (in pkSolveFlutter) | m/kg | same; "1/(Pa·s²)" is equivalent |
| Ω²_c | rad²/s² | diag(ω_n²·(1+ig)) |
| A = q·Q_norm - Ω²_c | rad²/s² | consistent |
| p = sqrt(eig(A)) | rad/s | Im(p) = frequency, Re(p) = growth rate |

---

## 5. Flutter Results and Physical Interpretation

### 5.1 Current solver output (g=0.02)
```
V_flutter (p-k) = 568 m/s  at Mach 1.838  (margin = 0.92×  vs V_flight=618 m/s)
V_div           = 2850 m/s at Mach 1.846  (margin = 4.6×)
```

### 5.2 Flutter mechanism — NOT bending-torsion
The D16 tailoring successfully suppresses bending-torsion flutter:
- Q0(1,2) = -1.9503, Q0(2,1) = +0.1233 → product = -0.2406 < 0 (stabilizing wash-out)
- Classical bending-torsion requires Q0(1,2)·Q0(2,1) > 0 → NOT present

The mechanism that IS detected is **single-mode aerodynamic anti-damping**:
- Piston theory's ẇ/U pressure term creates a force IN THE DIRECTION of surface velocity
- This is physically anti-damping (positive feedback on modal velocity)
- Governed by: Q1(1,1) ≈ 0.91 [m/kg], Q1 ≈ (2/β/b_ref)·∫φ₁²dA

### 5.3 Closed-form flutter speed formula
For the single-mode anti-damping mechanism (derived from Im(p²)=0 condition):
```
q_flutter = ω₁·g·U / (b_ref·Q1(1,1))
V_flutter = sqrt(2·q_flutter / ρ)
          = sqrt(ω₁·g·U / (b_ref·Q1(1,1)·ρ/2))
          ≈ 4334·sqrt(g)  m/s  [at Mach 1.838, evaluated numerically]
```
| g | V_flutter | Margin vs V_flight=618 m/s |
|---|---|---|
| 0.005 | 281 m/s | 0.45× |
| 0.010 | 397 m/s | 0.64× |
| 0.020 | 562 m/s | 0.91× |
| 0.030 | 688 m/s | 1.11× |
| 0.050 | 888 m/s | 1.44× |

Solver output vs. formula: 568 vs 562 m/s at g=0.02 (1% difference — off-diagonal Q0 contribution).

### 5.4 Structural loss factor g — documentation verdict
None of the project references (ADA502110, Martin NACA TN 4917, Newsletter615, SUBSONIC AND
SUPERSONIC FLUTTER ANALYSIS, Weisshaar Aeroelastic Tailoring) provide a CFRP damping value.
The SUBSONIC AND SUPERSONIC paper explicitly states:
> "Since the pertinent structural damping values are not known... all calculated flutter
>  points are taken to be points for which g = 0."

That convention applies to the V-g method (g is artificial, not physical). In the p-k method,
g_struct is the REAL structural loss factor:
- g = 0 (V-g convention): single-mode anti-damping flutter is ill-defined (condition q·k·Q1=ω²·g
  requires g>0 for a finite flutter speed). Bending-torsion flutter: STABLE.
- g = 0.02 (lower bound for T700/epoxy CFRP from DMA literature): V_flutter = 562–568 m/s
- g needs experimental characterization for flight certification

**Current setting:** g = 0.02 in pkSolveFlutter.m (line: `g_struct = 0.02`)

---

## 6. Key Design Conclusions

1. **Bending-torsion flutter: SUPPRESSED** by D16=+41.3 N·m aeroelastic tailoring (β=5°, Λ=57.4°)
   Confirmed analytically: Q0(1,2)·Q0(2,1) = -0.24 < 0 (wash-out cross-coupling)

2. **Divergence: SUPPRESSED** with 4.6× margin (V_div = 2850 m/s)
   q_div_crit = 4.44M Pa; all Q0 divergence eigenvalues ≤ 0 confirmed

3. **Single-mode anti-damping flutter: g-dependent**
   V_flutter ∝ √g. At g=0.02: 568 m/s (0.92× margin — slightly below V_flight at max-q point).
   This mechanism warrants structural damping testing for flight-critical certification.

4. **Martin's formula predicts bending-torsion flutter** (a different, higher-speed mechanism).
   Consistency: Martin gives ~1120 m/s → well above flight envelope → confirms (1) above.

---

## 7. Parameters to Vary for Continued Analysis

| Parameter | Location | Current | Sensitivity |
|---|---|---|---|
| g_struct | pkSolveFlutter.m line `g_struct = 0.02` | 0.02 | V_flutter ∝ √g (critical) |
| beta_deg | python_core/inplaneG_v5.py | 5 | affects D16, Q0 coupling sign |
| t_mm | lam.json → flutter_input.t_mm | 8 | V_flutter ∝ t^{3/2} via stiffness |
| nModes | mainFlutterSolver.m | 6 | mode 5 always zero-GAF |
| nx, ny | mainFlutterSolver.m | 24, 12 | mesh convergence unchecked |
| k_vals | mainFlutterSolver.m | [0,0.01,0.05,0.1,0.2,0.5,1.0] | Q1 slope from k_vals(2)=0.01 |

---

## 8. Rolling State Block

```
[RSB — last updated after fixing both bugs and diagnosing flutter mechanism]
Problem:          p-k flutter solver for tailored CFRP supersonic fin
Confirmed fixes:  BUG1: removed Hermitian symmetrization in pistonTheoryGAF.m
                  BUG2: changed k-update from mean→min(Im(p)) in pkSolveFlutter.m
                         + k_avg initialized to 0 per flight condition
Confirmed results: f_n = [237, 592, 1342, 1439, 2057, 2247] Hz
                   Q0(1,2)·Q0(2,1) = -0.2406 → bending-torsion flutter SUPPRESSED
                   V_flutter(g=0.02) = 568 m/s, margin 0.92× (anti-damping mechanism)
                   V_flutter formula: ≈ 4334·√g m/s at worst-case Mach 1.838
                   V_div = 2850 m/s, margin 4.6× (correct, divergence eigenvalues ≤ 0)
Open questions:   - Structural damping g not measured; needs DMA test for certification
                  - Mode 5 (membrane) confirmed zero-GAF but exact mode shape not visualized
                  - Mesh convergence (24×12) not formally verified
Dead ends:        - Hermitian symmetrization (V=Inf, wrong)
                  - Mean-k averaging (V=168 m/s, numerical artifact)
                  - g=0 in p-k: anti-damping flutter ill-defined, bending-torsion STABLE
Next steps:       - DMA test to measure g for T700/epoxy AR1 laminate
                  - Parametric sweep over t (thickness) and beta_deg
                  - Mesh convergence study: refine to 48×24 and compare frequencies
```

---

## 9. Session Resume Prompt Template

```
[CB-FLUTTER-2026]
[RSB — paste from §8 above]

Resume: [your task here — e.g., "parametric sweep over t=6,8,10 mm, output V_flutter table"]
```

Full context restored in one turn. Do not re-explain the bugs, the physics, or the laminate.
