"""
inplaneG_v5.py 
====================================================================
Classical Laminate Theory (CLPT) + Halpin-Tsai Micromechanics
Aeroelastic tailoring via global beta-rotation (Weisshaar 1981)

Upgrades over v4
----------------
  * CLI  --vf, --beta, --layup, --thickness, --json, --quiet
  * JSON export pipeline for flutterEstimate_v3.py
  * Target-thickness scaling (proportional t-scaling)
  * Four half-stack architectures: ar1, more_db, more_ga, equal
  * B-matrix symmetry check (relative tolerance)

References
----------
  Halpin & Tsai (1969)  AFML-TR-67-423
  Jones (1999)  Mechanics of Composite Materials  2nd ed.
  Weisshaar (1981)  J. Aircraft 18(8):669-676
  Shirk, Hertz & Weisshaar (1986)  J. Aircraft 23(1):6-21
  Naik & Shembekar (1992)  J. Compos. Mater. 26(15):2196

Usage
-----
  python3 inplaneG_v5.py                         # AR1 layup, full beta sweep
  python3 inplaneG_v5.py --beta 20               # single beta value
  python3 inplaneG_v5.py --vf 0.55 --beta 20
  python3 inplaneG_v5.py --thickness 6.0         # scale to 6 mm total
  python3 inplaneG_v5.py --json lam.json         # export for flutter script
  python3 inplaneG_v5.py --quiet --json lam.json # silent JSON export
"""

import math
import sys
import json
import argparse


# ============================================================
# CLI
# ============================================================
def parse_args():
    p = argparse.ArgumentParser(
        description="CLPT Laminate Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    p.add_argument("--vf",        type=float, default=0.50,
                   help="Fiber volume fraction (default 0.50)")
    p.add_argument("--beta",      type=float, default=None,
                   help="Single beta offset angle [deg]. Omit for full sweep.")
    p.add_argument("--layup",     choices=["ar1","more_db","more_ga","equal"],
                   default="ar1", help="Half-stack architecture (default ar1)")
    p.add_argument("--thickness", type=float, default=None,
                   help="Target total laminate thickness [mm]. Scales all plies.")
    p.add_argument("--json",      type=str,   default=None,
                   help="Write results to this JSON file (for flutter pipeline)")
    p.add_argument("--quiet",     action="store_true",
                   help="Suppress console output (useful with --json)")
    return p.parse_args()


# ============================================================
# MATERIAL DATABASE  —  T700 Carbon / Epoxy
# ============================================================
Ef1   = 230e9;  Ef2  = 15e9;   Gf12 = 27e9;  nuf  = 0.20
Em    = 3.5e9;  num  = 0.35;   Gm   = Em / (2.0 * (1.0 + num))
rho_f = 1800.0   # kg/m3

FAW_DB300 = 0.300   # kg/m2 — CARBONODB300 +-45 NCF biaxial
FAW_GA90R = 0.302   # kg/m2 — CARBONOGA90R woven 0/90
kc        = 0.92    # crimp knockdown for woven (Naik & Shembekar 1992)


def halpin_tsai(Ep, Em_, Vf, xi):
    eta = ((Ep / Em_) - 1.0) / ((Ep / Em_) + xi)
    return Em_ * (1.0 + xi * eta * Vf) / (1.0 - eta * Vf)


def material_props(Vf):
    """
    Return per-ply elastic constants and thicknesses for both fabrics.

    CARBONODB300 (NCF +-45): no crimp penalty; UD constants apply directly.
      Two independent plies per physical layer — each entered separately.
      t_single = (FAW/2) / (Vf * rho_f)

    CARBONOGA90R (plain woven 0/90): balanced => E1 = E2.
      One physical layer = one stack entry. Crimp knockdown kc = 0.92.
      t_ply = FAW / (Vf * rho_f)
    """
    E1   = Ef1 * Vf + Em * (1.0 - Vf)       # ROM — fiber-dominated
    nu12 = nuf * Vf  + num * (1.0 - Vf)      # ROM
    E2   = halpin_tsai(Ef2,  Em, Vf, xi=2.0)
    G12  = halpin_tsai(Gf12, Gm, Vf, xi=1.0)

    t_DB = (FAW_DB300 / 2.0) / (Vf * rho_f)
    t_GA = FAW_GA90R / (Vf * rho_f)

    E1_GA = 0.5 * (E1 + E2) * kc   # cross-ply average + crimp
    # nu12 for balanced woven [0/90]: CLT cross-ply average of UD plies
    # nu12_woven = 2*Q12 / (Q11+Q22) = 2*nu12_UD*E2 / (E1+E2)  [exact]
    # Derivation: Jones (1999) §2 CLT applied to [0/90] stack; exact (not approx), ≈ 0.032 at Vf=0.50
    nu12_GA = 2.0 * nu12 * E2 / (E1 + E2)  # [Jones, Mechanics of Composite Materials, 1999 §2]
    return {
        'ud': (E1, E2, G12, nu12),
        'DB': (E1,    E2, G12, nu12, t_DB),       # NCF: UD constants
        'GA': (E1_GA, E1_GA, G12 * kc, nu12_GA, t_GA)  # woven: E1=E2; kc on G12 — [Naik & Shembekar 1992]; nu12 — [Jones 1999]
    }


# ============================================================
# FIN GEOMETRY
# ============================================================
SWEEP_DEG = 57.4


# ============================================================
# CLPT CORE
# ============================================================
def Q_reduced(e1, e2, g12, nu12):
    nu21 = nu12 * e2 / e1
    D = 1.0 - nu12 * nu21
    return e1/D, e2/D, nu12*e2/D, g12


def Qbar(Q11, Q22, Q12, Q66, theta_deg):
    """Off-axis Qbar via standard CLPT 4th-order rotation."""
    a  = math.radians(theta_deg)
    m, n = math.cos(a), math.sin(a)
    m2, n2, mn = m*m, n*n, m*n
    return (
        Q11*m2*m2 + 2*(Q12+2*Q66)*n2*m2 + Q22*n2*n2,   # 11
        Q11*n2*n2 + 2*(Q12+2*Q66)*n2*m2 + Q22*m2*m2,   # 22
        (Q11+Q22-4*Q66)*n2*m2 + Q12*(m2*m2+n2*n2),      # 12
        (Q11+Q22-2*Q12-2*Q66)*n2*m2 + Q66*(n2*n2+m2*m2),# 66
        (Q11-Q12-2*Q66)*m2*mn - (Q22-Q12-2*Q66)*n2*mn,  # 16
        (Q11-Q12-2*Q66)*n2*mn - (Q22-Q12-2*Q66)*m2*mn,  # 26
    )


def build_ABD(stack):
    """
    Integrate A, B, D matrices through thickness.
    stack: list of (angle_deg, e1, e2, g12, nu12, t)
    Returns dict with all ABD terms + G_xy, G_eff, t_total.
    """
    t_total = sum(p[5] for p in stack)
    z = -t_total / 2.0
    A = [0.0]*6;  B = [0.0]*6;  D = [0.0]*6   # [11,22,12,66,16,26]
    for ang, e1, e2, g12, nu12, t in stack:
        q = Q_reduced(e1, e2, g12, nu12)
        qb = Qbar(*q, ang)
        z0, z1 = z, z + t
        dz  = z1 - z0
        dz2 = (z1**2 - z0**2) / 2.0
        dz3 = (z1**3 - z0**3) / 3.0
        for i, q_ in enumerate(qb):
            A[i] += q_ * dz
            B[i] += q_ * dz2
            D[i] += q_ * dz3
        z = z1
    G_xy  = A[3] / t_total
    G_eff = 12.0 * D[3] / t_total**3
    keys  = ('A11','A22','A12','A66','A16','A26',
             'B11','B22','B12','B66','B16','B26',
             'D11','D22','D12','D66','D16','D26')
    vals  = A + B + D
    return {**dict(zip(keys, vals)), 'G_xy': G_xy, 'G_eff': G_eff, 't_total': t_total}


# ============================================================
# STACK BUILDERS
# ============================================================
_FAB = {"B": "DB", "W": "GA"}

def make_ply(props, fabric, angle):
    e1, e2, g12, nu12, t = props[_FAB[fabric]]
    return (angle, e1, e2, g12, nu12, t)


def make_symmetric(half_desc, beta, props):
    """
    Symmetric balanced laminate from a half-stack descriptor.

    half_desc : list of ('B'|'W', base_angle_deg)
                NCF 'B' entries must appear as explicit ±45 pairs so the
                laminate is balanced (A16=A26=0). The bottom half is the
                exact mirror of the top (reversed order, same angles), giving
                B=0 by construction.
    beta      : global tailoring rotation [deg] — applied uniformly to all plies.

    Physics (Weisshaar 1981):
      beta=0  => D16≠0 from ±45 z-asymmetry within each pair
      beta=20 => |D16| amplified — washout for aft-swept fin (stabilising)
      beta=45 => D16=0 (NCF rotates to 0°/90°; woven rotation-neutral at E1=E2)
      beta>45 => D16 reverses sign — wash-in, destabilising
    """
    top = [(fk, ang + beta) for fk, ang in half_desc]
    bot = [(fk, ang + beta) for fk, ang in reversed(half_desc)]
    return [make_ply(props, fk, ang) for fk, ang in top + bot]


def scale_thickness(stack, target_m):
    """Scale all ply thicknesses proportionally to reach target_m total."""
    cur = sum(p[5] for p in stack)
    f   = target_m / cur if cur > 0 else 1.0
    return [(a, e1, e2, g12, nu, t*f) for a, e1, e2, g12, nu, t in stack]


# ============================================================
# HALF-STACK LIBRARY
# ============================================================
LAYUPS = {
    # AR1 All-Rounder: +-45 NCF outer skins + 0/90 woven backbone
    'ar1': [
        ('B',+45),('B',-45),
        ('W',  0),('W',  0),
        ('B',+45),('B',-45),
        ('W',  0),
        ('B',+45),('B',-45),
        ('W',  0),('W',  0),
    ],
    # Torsion-biased: more NCF for D66
    'more_db': [
        ('B',+45),('B',-45),('B',+45),('B',-45),
        ('W',  0),('W',  0),
        ('B',+45),('B',-45),
        ('W',  0),
        ('B',+45),('B',-45),
        ('W',  0),
        ('B',+45),('B',-45),
    ],
    # Bending-biased: more 0/90 woven for D11
    'more_ga': [
        ('W',  0),('W',  0),
        ('B',+45),('B',-45),
        ('W',  0),('W',  0),
        ('B',+45),('B',-45),
        ('W',  0),('W',  0),
        ('B',+45),('B',-45),
        ('W',  0),('W',  0),
    ],
    # Equal mix
    'equal': [
        ('B',+45),('B',-45),
        ('W',  0),('W',  0),
        ('B',+45),('B',-45),
        ('W',  0),('W',  0),
        ('B',+45),('B',-45),
        ('W',  0),('W',  0),
    ],
}


# ============================================================
# PRINTERS
# ============================================================
W = "=" * 70
D = "-" * 70


def print_material(Vf, props):
    E1, E2, G12, nu12 = props['ud']
    _, _, _, _, t_DB  = props['DB']
    E1g, _, _, _, t_GA = props['GA']
    print(W)
    print("  MATERIAL CARDS")
    print(W)
    print(f"\n  Vf = {Vf:.2f}")
    print(f"\n  UD T700/Epoxy (Halpin-Tsai micromechanics)")
    print(f"    E1   = {E1/1e9:.3f} GPa   E2  = {E2/1e9:.3f} GPa")
    print(f"    G12  = {G12/1e9:.3f} GPa   nu12 = {nu12:.4f}")
    print(f"\n  CARBONODB300  +-45 NCF biaxial  (FAW={FAW_DB300*1000:.0f} g/m2)")
    print(f"    t_single = {t_DB*1e3:.4f} mm/ply   t_pair = {t_DB*2e3:.4f} mm")
    print(f"    UD constants apply (no crimp)")
    print(f"\n  CARBONOGA90R  woven 0/90  (FAW={FAW_GA90R*1000:.0f} g/m2, 12K, kc={kc})")
    print(f"    t_ply  = {t_GA*1e3:.4f} mm")
    print(f"    E1=E2  = {E1g/1e9:.3f} GPa  (cross-ply avg + crimp knockdown)")
    print(f"    Note: E1=E2 => D11 invariant under 90 deg rotation for this fabric")


def print_beta_sweep(half_desc, beta_list, props, target_t):
    print(f"\n{W}")
    print(f"  AEROELASTIC TAILORING  —  beta SWEEP  (Weisshaar 1981)")
    print(W)
    print(f"  {'beta':>6}  {'D11':>9}  {'D22':>9}  {'D66':>9}  {'D16':>10}  {'G_xy':>9}  {'t':>6}  Note")
    print(f"  {'[deg]':>6}  {'[N.m]':>9}  {'[N.m]':>9}  {'[N.m]':>9}  {'[N.m]':>10}  {'[GPa]':>9}  {'[mm]':>6}")
    print("  " + D)
    results = []
    for b in beta_list:
        stk = make_symmetric(half_desc, b, props)
        if target_t:
            stk = scale_thickness(stk, target_t * 1e-3)
        r = build_ABD(stk)
        note = ""
        if   b == 0:                         note = "standard"
        elif b == 15:                        note = "tailoring onset"
        elif b == 20:                        note = "RECOMMENDED — washout"
        elif b == 30:                        note = "peak |D16|"
        elif b == 45:                        note = "D16=0 (symmetry zero)"
        elif abs(b - SWEEP_DEG) < 1.5:      note = f"aligned with sweep ({SWEEP_DEG} deg) — DESTABILISING"
        elif b == 90:                        note = "same D66 as beta=0; D16 sign reversed"
        print(f"  {b:>6.1f}  {r['D11']:>9.2f}  {r['D22']:>9.2f}  {r['D66']:>9.2f}  "
              f"{r['D16']:>10.4f}  {r['G_xy']/1e9:>9.4f}  {r['t_total']*1e3:>6.3f}  {note}")
        results.append({'beta': b, **r})
    print(f"""
  Stabilising range for aft-swept fin (Lambda={SWEEP_DEG} deg): 0 < beta < 45 deg
  At beta=20 deg: D16 < 0  =>  WASHOUT
    Mechanism: tip bends up => D16 torques TE up => AoA decreases
               => aero load drops => negative flutter feedback => stable
  At beta >= 45 deg (especially near 57 deg): D16 > 0  => WASH-IN  => DESTABILISING
  Reference: Weisshaar (1981), Shirk Hertz & Weisshaar (1986)
""")
    return results


def print_layup(label, stk, r):
    print(f"\n{W}")
    print(f"  LAYUP REPORT  —  {label}")
    print(W)
    print(f"\n  {'#':<4} {'Fabric':<18} {'Angle':>9}  {'t [mm]':>8}  {'z_mid [mm]':>11}")
    print("  " + D)
    z = -r['t_total'] / 2.0
    n_db = 0; n_ga = 0
    for i, (ang, e1, e2, g12, nu12, t) in enumerate(stk):
        is_db  = e1 > 100e9
        fabric = 'DB300 NCF +-45' if is_db else 'GA90R wov 0/90'
        z_mid  = z + t / 2.0
        sign   = '+' if ang >= 0 else ''
        print(f"  {i+1:<4} {fabric:<18}  {sign}{ang:+6.1f} deg  {t*1e3:>8.4f}  {z_mid*1e3:>11.4f}")
        z += t
        if is_db: n_db += 1
        else:     n_ga += 1
    print(f"\n  Ply count: {len(stk)} total  ({n_db} DB300 singles + {n_ga} GA90R layers)")
    print(f"\n  ABD MATRIX — KEY TERMS")
    print(f"    D11  = {r['D11']:>10.4f} N.m   spanwise bending stiffness")
    print(f"    D22  = {r['D22']:>10.4f} N.m   chordwise bending stiffness")
    print(f"    D66  = {r['D66']:>10.4f} N.m   torsional stiffness (flutter driver)")
    print(f"    D16  = {r['D16']:>10.4f} N.m   bend-twist coupling (tailoring)")
    print(f"    D12  = {r['D12']:>10.4f} N.m")
    print(f"    G_xy = {r['G_xy']/1e9:>10.4f} GPa  in-plane shear")
    print(f"    Geff = {r['G_eff']/1e9:>10.4f} GPa  12*D66/t^3  (flutter input)")
    print(f"    t    = {r['t_total']*1e3:>10.4f} mm")
    # Symmetry check — relative B norm
    b_max = max(abs(r[k]) for k in ('B11','B22','B12','B66','B16','B26'))
    # Reference: D66/t [N.m/m = N] matches B dimensions [N] — dimensionally appropriate.
    # Using A66*t overestimates reference for thick laminates. ERRORS.md I-3.
    b_ref = abs(r['D66']) / r['t_total']  # [N.m] / [m] = [N]
    b_rel = b_max / b_ref if b_ref > 0 else b_max
    ok    = b_rel < 1e-5
    print(f"\n  SYMMETRY CHECK (B=0 for symmetric laminate)")
    print(f"    max|B_ij| = {b_max:.2e} N   relative = {b_rel:.2e}  "
          f"{'OK' if ok else 'FAIL — check stacking'}")

# ============================================================
# MAIN
# ============================================================
def main():
    args       = parse_args()
    Vf         = args.vf
    target_t   = args.thickness          # mm or None
    props      = material_props(Vf)
    half_desc  = LAYUPS[args.layup]
    quiet      = args.quiet

    if not quiet:
        print_material(Vf, props)

    beta_list  = ([args.beta] if args.beta is not None
                  else [0, 15, 20, 30, 45, round(SWEEP_DEG), 75, 90])

    if not quiet:
        print_beta_sweep(half_desc, beta_list, props, target_t)

    # Detailed layup reports — always beta=0 baseline; tailored only when --beta is given
    beta_pairs = [("STANDARD  beta=0 deg", 0.0)]
    if args.beta is not None:
        beta_pairs.append((f"TAILORED  beta={args.beta:.1f} deg", args.beta))

    reports = {}
    for label, b in beta_pairs:
        stk = make_symmetric(half_desc, b, props)
        if target_t:
            stk = scale_thickness(stk, target_t * 1e-3)
        r = build_ABD(stk)
        reports[b] = (stk, r)
        if not quiet:
            print_layup(label, stk, r)

    r_std  = reports[0.0][1]
    r_tail = reports[args.beta][1] if args.beta is not None else r_std

    # JSON export
    if args.json:
        E1, E2, G12, nu12 = props['ud']
        _, _, _, _, t_DB  = props['DB']
        E1g, _, _, _, t_GA = props['GA']
        out = {
            "tool": "inplaneG_v5",
            "inputs": {"Vf": Vf, "layup": args.layup, "target_thickness_mm": target_t},
            "material": {
                "fiber": "T700 12K", "matrix": "Epoxy LY1564",
                "Vf": Vf,
                "E1_ud_GPa": E1/1e9, "E2_ud_GPa": E2/1e9,
                "G12_ud_GPa": G12/1e9, "nu12_ud": nu12,
                "t_DB_single_mm": t_DB*1e3, "t_GA_mm": t_GA*1e3,
            },
            "standard_beta0": {
                "D11_Nm": r_std['D11'], "D22_Nm": r_std['D22'],
                "D66_Nm": r_std['D66'], "D16_Nm": r_std['D16'],
                "G_xy_GPa": r_std['G_xy']/1e9,
                "G_eff_GPa": r_std['G_eff']/1e9,
                "t_total_mm": r_std['t_total']*1e3,
            },
            f"tailored_beta": {
                "D11_Nm": r_tail['D11'], "D22_Nm": r_tail['D22'],
                "D66_Nm": r_tail['D66'], "D16_Nm": r_tail['D16'],
                "G_xy_GPa": r_tail['G_xy']/1e9,
                "G_eff_GPa": r_tail['G_eff']/1e9,
                "t_total_mm": r_tail['t_total']*1e3,
            },
            "flutter_input": {
                "D66_conservative_Nm": r_tail['D66'],
                "t_mm": r_tail['t_total']*1e3,
                "note": "Use tailored D66 (beta=20) — conservative lower bound"
            }
        }
        with open(args.json, 'w') as f:
            json.dump(out, f, indent=2)
        if not quiet:
            print(f"\n  JSON exported -> {args.json}")

    return r_std, r_tail


if __name__ == "__main__":
    main()
