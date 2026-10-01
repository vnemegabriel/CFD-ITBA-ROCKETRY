"""
inplaneG_v5.py 
====================================================================
Classical Laminate Theory (CLPT) + Halpin-Tsai Micromechanics
Aeroelastic tailoring via global beta-rotation (Weisshaar 1981)

"""

import math
import sys
import json
import argparse
import itertools


# ============================================================
# CLI
# ============================================================
def parse_args():
    p = argparse.ArgumentParser(
        description="CLPT Laminate Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    p.add_argument("--vf",        type=float, default=0.50,
                   help="Nominal fiber volume fraction (default: 0.50)")
    p.add_argument("--beta",      type=float, default=5,
                   help="Single beta offset angle [deg]. Omit for full sweep (default: 5)")
    p.add_argument("--layup",     choices=["ar1","more_db","more_ga","equal"],
                   default="ar1", help="Half-stack architecture (default: ar1)")
    p.add_argument("--thickness", type=float, default=8,
                   help="Target mold thickness [mm]. Code calculates required plies and true Vf")
    p.add_argument("--json",       type=str,   default=None,
                   help="Write results to this JSON file (default: lam<t>mm.json)")
    p.add_argument("--quiet",      action="store_true", default=False,
                   help="Suppress console output (useful with --json) (default: False)")
    p.add_argument("--sweep",      action="store_true", default=False,
                   help="Sweep beta 0 to 45 deg in steps of 5, export lam_sweep.json")
    p.add_argument("--sweep-json", type=str,   default="lam_sweep.json",
                   help="Output path for sweep JSON (default: lam_sweep.json)")
    args = p.parse_args()
    if args.json is None:
        t_str = str(args.thickness).replace(".", "_")
        args.json = f"lam{t_str}mm.json"
    return args


# ============================================================
# MATERIAL DATABASE  —  T700 Carbon / Epoxy
# ============================================================
Ef1   = 230e9;  Ef2  = 15e9;   Gf12 = 27e9;  nuf  = 0.20
Em    = 3.5e9;  num  = 0.35;   Gm   = Em / (2.0 * (1.0 + num))
rho_f = 1800.0   # kg/m³ — T700 carbon fiber
rho_m = 1200.0   # kg/m³ — cured epoxy (LY1564 / Huntsman)

FAW_DB300 = 0.300   # kg/m2 — CARBONODB300 +-45 NCF biaxial
FAW_GA90R = 0.302   # kg/m2 — CARBONOGA90R woven 0
kc        = 0.92    # crimp knockdown for woven (Naik & Shembekar 1992)


def halpin_tsai(Ep, Em_, Vf, xi):
    eta = ((Ep / Em_) - 1.0) / ((Ep / Em_) + xi)
    return Em_ * (1.0 + xi * eta * Vf) / (1.0 - eta * Vf)


def material_props(Vf):
    E1   = Ef1 * Vf + Em * (1.0 - Vf)       # ROM — fiber-dominated
    nu12 = nuf * Vf  + num * (1.0 - Vf)     # ROM
    E2   = halpin_tsai(Ef2,  Em, Vf, xi=2.0)
    G12  = halpin_tsai(Gf12, Gm, Vf, xi=1.0)

    t_DB = FAW_DB300 / 2.0 / (Vf * rho_f)
    t_GA = FAW_GA90R / (Vf * rho_f)

    # Composite density: rule of mixtures  ρ_c = ρ_f·Vf + ρ_m·(1-Vf)
    rho_c = rho_f * Vf + rho_m * (1.0 - Vf)

    E1_GA = 0.5 * (E1 + E2) * kc   # cross-ply average + crimp
    nu12_GA = 2.0 * nu12 * E2 / (E1 + E2) 
    
    # Rule of mixtures: ρ_lam = Vf·ρ_f + (1−Vf)·ρ_m  [uniform across both fabric types]
    rho_lam = Vf * rho_f + (1.0 - Vf) * rho_m

    return {
        'ud':      (E1, E2, G12, nu12),
        'DB':      (E1,    E2, G12, nu12, t_DB),
        'GA':      (E1_GA, E1_GA, G12 * kc, nu12_GA, t_GA),
        'rho_lam': rho_lam,   # kg/m³ — effective laminate density
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
    a  = math.radians(theta_deg)
    m, n = math.cos(a), math.sin(a)
    m2, n2, mn = m*m, n*n, m*n
    return (
        Q11*m2*m2 + 2*(Q12+2*Q66)*n2*m2 + Q22*n2*n2,   
        Q11*n2*n2 + 2*(Q12+2*Q66)*n2*m2 + Q22*m2*m2,   
        (Q11+Q22-4*Q66)*n2*m2 + Q12*(m2*m2+n2*n2),     
        (Q11+Q22-2*Q12-2*Q66)*n2*m2 + Q66*(n2*n2+m2*m2),
        (Q11-Q12-2*Q66)*m2*mn - (Q22-Q12-2*Q66)*n2*mn,  
        (Q11-Q12-2*Q66)*n2*mn - (Q22-Q12-2*Q66)*m2*mn,  
    )


def build_ABD(stack):
    t_total = sum(p[5] for p in stack)
    z = -t_total / 2.0
    A = [0.0]*6;  B = [0.0]*6;  D = [0.0]*6   
    for p in stack:
        ang, e1, e2, g12, nu12, t = p[:6]
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
    return (angle, e1, e2, g12, nu12, t, _FAB[fabric])


def make_symmetric(half_desc, beta, props):
    top = [(fk, ang + beta) for fk, ang in half_desc]
    bot = [(fk, ang + beta) for fk, ang in reversed(half_desc)]
    return [make_ply(props, fk, ang) for fk, ang in top + bot]


# ============================================================
# HALF-STACK LIBRARY
# ============================================================
LAYUPS = {
    'ar1': [
        ('B',+45),('B',-45),
        ('W',  0),('W',  0),
        ('B',+45),('B',-45),
        ('W',  0),
        ('B',+45),('B',-45),
        ('W',  0),('W',  0),
    ],
    'more_db': [
        ('B',+45),('B',-45),('B',+45),('B',-45),
        ('W',  0),('W',  0),
        ('B',+45),('B',-45),
        ('W',  0),
        ('B',+45),('B',-45),
        ('W',  0),
        ('B',+45),('B',-45),
    ],
    'more_ga': [
        ('W',  0),('W',  0),
        ('B',+45),('B',-45),
        ('W',  0),('W',  0),
        ('B',+45),('B',-45),
        ('W',  0),('W',  0),
        ('B',+45),('B',-45),
        ('W',  0),('W',  0),
    ],
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
    _, _, _, _, t_DB   = props['DB']
    E1g, _, _, _, t_GA = props['GA']
    rho_lam = props['rho_lam']
    print(W)
    print("  MATERIAL CARDS")
    print(W)
    print(f"\n  Vf = {Vf:.3f}")
    print(f"\n  UD T700/Epoxy (Halpin-Tsai micromechanics)")
    print(f"    E1   = {E1/1e9:.3f} GPa   E2  = {E2/1e9:.3f} GPa")
    print(f"    G12  = {G12/1e9:.3f} GPa   nu12 = {nu12:.4f}")
    print(f"\n  CARBONODB300  +-45 NCF biaxial  (FAW={FAW_DB300*1000:.0f} g/m2)")
    print(f"    t_single = {t_DB*1e3:.4f} mm/ply   t_pair = {t_DB*2e3:.4f} mm")
    print(f"    UD constants apply (no crimp)")
    print(f"\n  CARBONOGA90R  woven 0  (FAW={FAW_GA90R*1000:.0f} g/m2, 12K, kc={kc})")
    print(f"    t_ply  = {t_GA*1e3:.4f} mm")
    print(f"    E1=E2  = {E1g/1e9:.3f} GPa  (cross-ply avg + crimp knockdown)")
    print(f"\n  LAMINATE DENSITY  (rule of mixtures: rho = Vf*rho_f + (1-Vf)*rho_m)")
    print(f"    rho_fiber  = {rho_f:.0f} kg/m3   rho_matrix = {rho_m:.0f} kg/m3")
    print(f"    rho_lam    = {rho_lam:.1f} kg/m3   (at Vf = {Vf:.3f})")


def print_beta_sweep(half_desc, beta_list, props):
    print(f"\n{W}")
    print(f"  AEROELASTIC TAILORING  —  beta SWEEP  (Weisshaar 1981)")
    print(W)
    print(f"  {'beta':>6}  {'D11':>9}  {'D22':>9}  {'D66':>9}  {'D16':>10}  {'G_xy':>9}  {'t':>6}  Note")
    print(f"  {'[deg]':>6}  {'[N.m]':>9}  {'[N.m]':>9}  {'[N.m]':>9}  {'[N.m]':>10}  {'[GPa]':>9}  {'[mm]':>6}")
    print("  " + D)
    results = []
    for b in beta_list:
        stk = make_symmetric(half_desc, b, props)
        r = build_ABD(stk)
        note = ""
        if   b == 0:                         note = "standard"
        elif b == 15:                        note = "tailoring onset"
        elif b == 20:                        note = "RECOMMENDED — washout"
        elif b == 30:                        note = "peak |D16|"
        elif b == 45:                        note = "D16=0 (symmetry zero)"
        elif abs(b - SWEEP_DEG) < 1.5:       note = f"aligned with sweep ({SWEEP_DEG} deg) — DESTABILISING"
        elif b == 90:                        note = "same D66 as beta=0; D16 sign reversed"
        print(f"  {b:>6.1f}  {r['D11']:>9.2f}  {r['D22']:>9.2f}  {r['D66']:>9.2f}  "
              f"{r['D16']:>10.4f}  {r['G_xy']/1e9:>9.4f}  {r['t_total']*1e3:>6.3f}  {note}")
        results.append({'beta': b, **r})
    return results


def print_layup(label, stk, r, rho_lam=None):
    print(f"\n{W}")
    print(f"  LAYUP REPORT  —  {label}")
    print(W)
    print(f"\n  {'#':<4} {'Fabric':<18} {'Angle [deg]':>16}  {'t [mm]':>8}  {'z_mid [mm]':>11}")
    print("  " + D)
    z = -r['t_total'] / 2.0
    n_db = 0; n_ga = 0
    
    i = 0
    ply_num = 1
    while i < len(stk):
        ang, e1, e2, g12, nu12, t, fab_type = stk[i]
        
        # Group NCF DB300 pairs for printing
        if fab_type == 'DB' and i + 1 < len(stk):
            ang2, e1_2, e2_2, g12_2, nu12_2, t2, fab_type2 = stk[i+1]
            if fab_type2 == 'DB':
                fabric = 'DB300 NCF biaxial'
                z_mid = z + (t + t2) / 2.0
                ang_str = f"[{ang:+g}, {ang2:+g}]"
                print(f"  {ply_num:<4} {fabric:<18}  {ang_str:>16}  {(t+t2)*1e3:>8.4f}  {z_mid*1e3:>11.4f}")
                z += (t + t2)
                n_db += 1
                ply_num += 1
                i += 2
                continue
        
        # Woven or fallback isolated layer
        fabric = 'DB300 NCF single' if fab_type == 'DB' else 'GA90R woven 0'
        z_mid  = z + t / 2.0
        ang_str = f"[{ang:+g}]"
        print(f"  {ply_num:<4} {fabric:<18}  {ang_str:>16}  {t*1e3:>8.4f}  {z_mid*1e3:>11.4f}")
        z += t
        if fab_type == 'DB': n_db += 1
        else:                n_ga += 1
        ply_num += 1
        i += 1

    print(f"\n  Ply count: {ply_num - 1} physical layers  ({n_db} DB300 pairs + {n_ga} GA90R layers)")
    print(f"\n  ABD MATRIX — KEY TERMS")
    print(f"    D11  = {r['D11']:>10.4f} N.m   spanwise bending stiffness")
    print(f"    D22  = {r['D22']:>10.4f} N.m   chordwise bending stiffness")
    print(f"    D66  = {r['D66']:>10.4f} N.m   torsional stiffness (flutter driver)")
    print(f"    D16  = {r['D16']:>10.4f} N.m   bend-twist coupling (tailoring)")
    print(f"    D26  = {r['D26']:>10.4f} N.m   twist-shear coupling  *** exported to JSON ***")
    print(f"    D12  = {r['D12']:>10.4f} N.m")
    print(f"    G_xy = {r['G_xy']/1e9:>10.4f} GPa  in-plane shear")
    print(f"    Geff = {r['G_eff']/1e9:>10.4f} GPa  12*D66/t^3  (flutter input)")
    print(f"    t    = {r['t_total']*1e3:>10.4f} mm")
    if rho_lam is not None:
        areal_mass = rho_lam * r['t_total']   # kg/m²
        print(f"    rho  = {rho_lam:>10.1f} kg/m3 (Vf*rho_f + (1-Vf)*rho_m)")
        print(f"    m/A  = {areal_mass*1e3:>10.4f} g/m2  (surface mass = rho*t)")

    b_max = max(abs(r[k]) for k in ('B11','B22','B12','B66','B16','B26'))
    b_ref = abs(r['D66']) / r['t_total']
    b_rel = b_max / b_ref if b_ref > 0 else b_max
    ok    = b_rel < 1e-5
    print(f"\n  SYMMETRY CHECK (B=0 for symmetric laminate)")
    print(f"    max|B_ij| = {b_max:.2e} N   relative = {b_rel:.2e}  "
          f"{'OK' if ok else 'FAIL — check stacking'}")

    A16, A26 = r['A16'], r['A26']
    balanced = (abs(A16) + abs(A26)) / (r['A11'] + 1e-30) < 1e-4
    print(f"\n  BALANCE CHECK (A16=A26=0 for in-plane balanced laminate)")
    print(f"    A16 = {A16/1e6:>10.4f} MPa*mm   A26 = {A26/1e6:>10.4f} MPa*mm  "
          f"{'OK' if balanced else 'UNBALANCED -- beta rotation breaks balance; A16!=0 causes extension-shear coupling'}")


# ============================================================
# MAIN
# ============================================================
def main():
    args       = parse_args()
    Vf         = args.vf
    props      = material_props(Vf)
    half_desc  = LAYUPS[args.layup]
    quiet      = args.quiet

    # --------------------------------------------------------
    # SIZING / THICKNESS TARGETING (CLOSED-MOLD PHYSICS)
    # --------------------------------------------------------
    if args.thickness is not None:
        target_m = args.thickness / 1000.0
        target_half_m = target_m / 2.0
        
        # 1. Group into physical plies to prevent splitting DB300 pairs
        physical_plies = []
        i = 0
        while i < len(half_desc):
            fk, ang = half_desc[i]
            if fk == 'B' and i+1 < len(half_desc) and half_desc[i+1][0] == 'B':
                physical_plies.append([half_desc[i], half_desc[i+1]])
                i += 2
            else:
                physical_plies.append([half_desc[i]])
                i += 1
                
        def get_thick(phys_ply):
            return sum(props[_FAB[f]][4] for f, a in phys_ply)
            
        new_half = []
        current_half_t = 0.0
        n_db_half = 0   # DB300 biaxial plies (pairs) per half
        n_ga_half = 0   # GA90R woven plies per half

        # 2. Cycle pattern using NOMINAL Vf to find discrete ply count
        for ply in itertools.cycle(physical_plies):
            pt = get_thick(ply)
            err_before = abs(target_half_m - current_half_t)
            err_after  = abs(target_half_m - (current_half_t + pt))
            is_db = (ply[0][0] == 'B')

            if current_half_t + pt > target_half_m:
                if err_after < err_before:
                    new_half.extend(ply)
                    current_half_t += pt
                    if is_db: n_db_half += 1
                    else:     n_ga_half += 1
                break
            else:
                new_half.extend(ply)
                current_half_t += pt
                if is_db: n_db_half += 1
                else:     n_ga_half += 1

        if not new_half:  # Failsafe for unusually thin targets
            new_half.extend(physical_plies[0])
            if physical_plies[0][0][0] == 'B': n_db_half += 1
            else:                               n_ga_half += 1

        half_desc = new_half

        # 3. CLOSED MOLD PHYSICS - Recalculate TRUE Vf based on discrete plies
        total_FAW_half = 0.0
        for fk, ang in half_desc:
            if fk == 'B': total_FAW_half += FAW_DB300 / 2.0
            elif fk == 'W': total_FAW_half += FAW_GA90R
            
        total_FAW_symmetric = total_FAW_half * 2.0
        actual_Vf = total_FAW_symmetric / (rho_f * target_m)

        if not quiet:
            n_phys_half = n_db_half + n_ga_half
            print(f"\n{W}")
            print(f"  CLOSED-MOLD THICKNESS SIZING")
            print(W)
            print(f"  Target mold thickness  = {args.thickness:.2f} mm")
            print(f"  Nominal target Vf      = {args.vf:.3f}")
            print(f"  Physical plies per half= {n_phys_half}  ({n_db_half} DB300 biaxial + {n_ga_half} GA90R woven)")
            print(f"  Total physical plies   = {n_phys_half*2} (symmetric)  [{n_db_half*2} DB300 + {n_ga_half*2} GA90R]")
            print(f"  Total Areal Weight     = {total_FAW_symmetric:.3f} kg/m2")
            print(f"  --------------------------------------------------")
            print(f"  ACTUAL MOLDED Vf       = {actual_Vf:.3f}")
            print(f"  --------------------------------------------------")

            if actual_Vf > 0.65:
                print("  >>> WARNING: Vf exceeds 65%. High probability of dry spots. Increase target thickness or remove plies.")
            elif actual_Vf < 0.40:
                print("  >>> WARNING: Vf is below 40%. Highly resin rich. Reduce target thickness or add plies.")
        
        # 4. OVERWRITE MATERIAL PROPERTIES WITH TRUE Vf
        Vf = actual_Vf
        props = material_props(Vf)

    if not quiet:
        # Prints with the finalized (Actual Molded) properties
        print_material(Vf, props)

    beta_list  = ([args.beta] if args.beta is not None
                  else [0, 15, 20, 30, 45, round(SWEEP_DEG), 75, 90])

    if not quiet:
        print_beta_sweep(half_desc, beta_list, props)

    beta_pairs = [("STANDARD  beta=0 deg", 0.0)]
    if args.beta is not None:
        beta_pairs.append((f"TAILORED  beta={args.beta:.1f} deg", args.beta))

    reports = {}
    for label, b in beta_pairs:
        stk = make_symmetric(half_desc, b, props)
        r = build_ABD(stk)
        reports[b] = (stk, r)
        if not quiet:
            print_layup(label, stk, r, rho_lam=props['rho_lam'])

    r_std  = reports[0.0][1]
    r_tail = reports[args.beta][1] if args.beta is not None else r_std

    # JSON export
    if args.json:
        E1, E2, G12, nu12  = props['ud']
        _, _, _, _, t_DB   = props['DB']
        E1g, _, _, _, t_GA = props['GA']

        def serialize_stack(stk):
            """Convert full symmetric stack to ply-by-ply list for JSON."""
            t_total = sum(p[5] for p in stk)
            z = -t_total / 2.0
            plies = []
            i = 0
            ply_num = 1
            while i < len(stk):
                ang, e1, e2, g12, nu12_, t, fab_type = stk[i]
                if fab_type == 'DB' and i + 1 < len(stk) and stk[i+1][6] == 'DB':
                    ang2, _, _, _, _, t2, _ = stk[i+1]
                    z_mid = z + (t + t2) / 2.0
                    plies.append({
                        "ply": ply_num,
                        "fabric": "DB300_NCF_biaxial",
                        "fiber_angles_from_chord_deg": [round(ang, 4), round(ang2, 4)],
                        "t_mm": round((t + t2) * 1e3, 5),
                        "z_mid_mm": round(z_mid * 1e3, 5),
                    })
                    z += t + t2
                    i += 2
                else:
                    z_mid = z + t / 2.0
                    fab_label = "DB300_NCF_single" if fab_type == 'DB' else "GA90R_woven_0_90"
                    plies.append({
                        "ply": ply_num,
                        "fabric": fab_label,
                        "fiber_angles_from_chord_deg": [round(ang, 4)],
                        "t_mm": round(t * 1e3, 5),
                        "z_mid_mm": round(z_mid * 1e3, 5),
                    })
                    z += t
                    i += 1
                ply_num += 1
            return plies

        stk_tail = reports[args.beta][0] if args.beta is not None else reports[0.0][0]

        out = {
            "tool": "inplaneG_v5",
            "inputs": {
                "Vf_nominal": args.vf,
                "Vf_actual_molded": Vf,
                "layup": args.layup,
                "target_thickness_mm": args.thickness,
                "beta_deg": args.beta,
            },
            "material": {
                "fiber": "T700 12K", "matrix": "Epoxy LY1564",
                "Vf": Vf,
                "E1_ud_GPa": E1/1e9, "E2_ud_GPa": E2/1e9,
                "G12_ud_GPa": G12/1e9, "nu12_ud": nu12,
                "t_DB_single_mm": t_DB*1e3, "t_GA_mm": t_GA*1e3,
                "rho_fiber_kgm3": rho_f, "rho_matrix_kgm3": rho_m,
                "rho_lam_kgm3": round(props['rho_lam'], 2),
            },
            "standard_beta0": {
                "D11_Nm": r_std['D11'], "D22_Nm": r_std['D22'],
                "D12_Nm": r_std['D12'], "D66_Nm": r_std['D66'],
                "D16_Nm": r_std['D16'], "D26_Nm": r_std['D26'],
                "G_xy_GPa": r_std['G_xy']/1e9,
                "G_eff_GPa": r_std['G_eff']/1e9,
                "t_total_mm": r_std['t_total']*1e3,
            },
            "tailored_beta": {
                "D11_Nm": r_tail['D11'], "D22_Nm": r_tail['D22'],
                "D12_Nm": r_tail['D12'], "D66_Nm": r_tail['D66'],
                "D16_Nm": r_tail['D16'], "D26_Nm": r_tail['D26'],
                "G_xy_GPa": r_tail['G_xy']/1e9,
                "G_eff_GPa": r_tail['G_eff']/1e9,
                "t_total_mm": r_tail['t_total']*1e3,
            },
            "flutter_input": {
                "D66_conservative_Nm": r_tail['D66'],
                "t_mm": args.thickness if args.thickness is not None else r_tail['t_total']*1e3,
                "rho_mat_kgm3": round(props['rho_lam'], 2),
                "note": "Use tailored D66 (beta=20) -- conservative lower bound; t_mm is mold dimension"
            },
            "layup_stacking": {
                "angle_reference": "fiber angle from chordwise (LE->TE = 0 deg), CCW positive",
                "symmetry": "symmetric about mid-plane (B=0)",
                "total_physical_plies": len(serialize_stack(stk_tail)),
                "plies": serialize_stack(stk_tail),
            },
        }
        with open(args.json, 'w') as f:
            json.dump(out, f, indent=2)
        if not quiet:
            print(f"\n  JSON exported -> {args.json}")

    if args.sweep:
        _run_beta_sweep(half_desc, props, Vf, args)

    return r_std, r_tail


# ============================================================
# BETA SWEEP  — structural pre-processor for betaSweepSolver.m
# ============================================================
def _run_beta_sweep(half_desc, props, Vf, args):
    """Sweep beta 0:5:45, print verification table, export lam_sweep.json."""
    BETAS = list(range(0, 46, 5))

    SEP = "=" * 70
    print(f"\n{SEP}")
    print(f"  BETA SWEEP  (structural pre-processing for betaSweepSolver.m)")
    print(SEP)
    print(f"  {'beta':>5}  {'D11':>8}  {'D22':>8}  {'D66':>8}  "
          f"{'D16':>9}  {'D26':>9}  {'k_weiss':>9}")
    print(f"  {'[deg]':>5}  {'[N.m]':>8}  {'[N.m]':>8}  {'[N.m]':>8}  "
          f"{'[N.m]':>9}  {'[N.m]':>9}  {'[-]':>9}")
    print("  " + "-" * 66)

    configs = []
    for b in BETAS:
        stk = make_symmetric(half_desc, b, props)
        r   = build_ABD(stk)
        k_w = r['D16'] / math.sqrt(r['D11'] * r['D66'])
        print(f"  {b:>5}  {r['D11']:>8.3f}  {r['D22']:>8.3f}  {r['D66']:>8.3f}  "
              f"{r['D16']:>9.4f}  {r['D26']:>9.4f}  {k_w:>9.5f}")
        configs.append({
            "beta_deg":  b,
            "k_weiss":   round(k_w, 7),
            "D11_Nm":    r['D11'],  "D22_Nm": r['D22'],
            "D12_Nm":    r['D12'],  "D66_Nm": r['D66'],
            "D16_Nm":    r['D16'],  "D26_Nm": r['D26'],
        })

    out = {
        "tool": "inplaneG_v5_sweep",
        "sweep_info": {
            "layup":            args.layup,
            "thickness_mm":     args.thickness,
            "Vf_actual_molded": Vf,
            "t_mm":             args.thickness,
            "rho_mat_kgm3":     round(props['rho_lam'], 2),
            "angle_reference":  "fiber angle from chordwise (LE->TE = 0 deg), CCW positive",
            "betas_deg":        BETAS,
        },
        "configs": configs,
    }

    path = getattr(args, 'sweep_json', 'lam_sweep.json')
    with open(path, 'w') as f:
        json.dump(out, f, indent=2)

    print(f"\n  Sweep JSON -> {path}")


if __name__ == "__main__":
    main()