"""
flutterEstimate_v4.py  —  Fin Flutter Boundary Analysis (Supersonic)
========================================================================
Theory  : Ackeret aerodynamics (M > 1) throughout.
          NACA TN 4197 K-factor for structural geometry scaling only.
          Sweep correction: Bisplinghoff, Ashley & Halfman (1955) §5.5

Single-equation derivation
--------------------------
Moment arms about the elastic axis x_ea:
    e_sub = x_ea - 0.25   (subsonic  AC at 0.25c)
    e_sup = 0.50 - x_ea   (supersonic AC at 0.50c, Ackeret)
    R     = π·e_sub / (2·e_sup)              [≈ 2.356 for x_ea = 0.40c]

Structural speed scale (NACA TN 4197 K-factor, geometry only):
    S² = G_eff · t³ · (tc_ref/tc) / (ρ · b² · K)

Self-consistent supersonic flutter condition:
    V_f² = S² · R · √(V_f²/a² − 1)

Substituting u = Mf² and squaring:
    u² − A²·u + A² = 0    where  A = (S/a)² · R

Physical root (requires A > 2):
    u⁺ = (A²/2) · (1 + √(1 − 4/A²))
    V_f,super = a · √u⁺
    V_f,swept = V_f,super / cos(Λ)

Usage
-----
  python3 flutterEstimate_v4.py
  python3 flutterEstimate_v4.py --json laminate.json
  python3 flutterEstimate_v4.py --d66 206.18 --t 5.356 --mach 1.942
  python3 flutterEstimate_v4.py --sweep-table

References
----------
[1] Martin, H.C.        NACA TN 4197  (1958)
[2] Ackeret, J.         NACA TM 317   (1925)
[3] Bisplinghoff, Ashley & Halfman.  Aeroelasticity (1955) §5.5
"""
import math, json, argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / "core" / "python"))
from fin_flutter.flight_data import read_flight_data

# ─── CLI ────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="Fin Flutter Estimation v4")
    p.add_argument("--json",        type=str,   default=None)
    p.add_argument("--d66",         type=float, default=None,  help="D66 [N.m]")
    p.add_argument("--t",           type=float, default=None,  help="thickness [mm]")
    p.add_argument("--altitude",    type=float, default=1462.0,help="max-q altitude [m]")
    p.add_argument("--velocity",    type=float, default=None,  help="rocket speed [m/s]")
    p.add_argument("--mach",        type=float, default=1.942, help="rocket Mach")
    p.add_argument("--span",        type=float, default=0.160, help="semi-span [m]")
    p.add_argument("--cr",          type=float, default=0.300, help="root chord [m]")
    p.add_argument("--ct",          type=float, default=0.150, help="tip chord [m]")
    p.add_argument("--sweep",       type=float, default=57.4,  help="LE sweep [deg]")
    p.add_argument("--x-ea",        type=float, default=0.40,  help="EA position [x/c]")
    p.add_argument("--sweep-off",   action="store_true")
    p.add_argument("--fsf-req",     type=float, default=1.50)
    p.add_argument("--sweep-table", action="store_true")
    p.add_argument("--csv",         type=str,   default=None)
    return p.parse_args()

# ─── ISA ATMOSPHERE ─────────────────────────────────────────────────────────
def isa(h):
    assert 0 <= h <= 11000, f"ISA troposphere only valid 0–11 km, got {h} m"
    T0, p0, L, g, R, gam = 288.15, 101325.0, 0.0065, 9.80665, 287.058, 1.4
    T   = T0 - L * h
    p   = p0 * (T / T0) ** (g / (R * L))
    rho = p / (R * T)
    a   = math.sqrt(gam * R * T)
    return rho, p, T, a

# ─── FIN GEOMETRY ───────────────────────────────────────────────────────────
def fin_ar(b, cr, ct):
    return b**2 / (0.5 * (cr + ct) * b)

def fin_mac(b, cr, ct):
    lam = ct / cr
    return cr * (2/3) * (1 + lam + lam**2) / (1 + lam)

def sweep_factor(sweep_deg):
    if sweep_deg == 0.0:  return 1.0
    if sweep_deg  > 0.0:  return 1.0 / math.cos(math.radians(sweep_deg))
    return math.cos(math.radians(abs(sweep_deg)))

# ─── SUPERSONIC QUADRATIC FLUTTER ───────────────────────────────────────────
def flutter_supersonic_quadratic(G_eff, rho, t, b, AR, tc, a,
                                 x_ea=0.40, sweep_deg=0.0):
    K      = max(0.65 * (1.0 + 0.10 * (AR - 1.0)), 0.40)
    S2     = G_eff * t**3 * (0.02 / tc) / (rho * b**2 * K)
    S      = math.sqrt(max(S2, 0.0))

    e_sub  = x_ea - 0.25
    e_sup  = 0.50  - x_ea
    if e_sub <= 0.0:
        raise ValueError(f"x_ea = {x_ea:.3f} ≤ 0.25: EA at or forward of subsonic AC.")
    if e_sup <= 0.0:
        raise ValueError(f"x_ea = {x_ea:.3f} ≥ 0.50: EA at or aft of supersonic AC.")

    R      = math.pi * e_sub / (2.0 * e_sup)
    A      = (S / a)**2 * R

    if A <= 2.0:
        raise ValueError(
            f"A = {A:.4f} ≤ 2.0: no real supersonic flutter root "
            f"(S={S:.1f} m/s, a={a:.2f} m/s, R={R:.4f}). Increase G_eff."
        )

    disc       = 1.0 - 4.0 / A**2
    u_plus     = (A**2 / 2.0) * (1.0 + math.sqrt(disc))
    Vf_super   = a * math.sqrt(u_plus)
    f_sw       = sweep_factor(sweep_deg)
    Vf_swept   = Vf_super * f_sw

    return dict(S=S, R=R, A=A, disc=disc, u_plus=u_plus,
                Vf_super=Vf_super, Mf_super=Vf_super/a,
                f_sw=f_sw, Vf_swept=Vf_swept, Mf_swept=Vf_swept/a)

# ─── COMPUTE ────────────────────────────────────────────────────────────────
def compute_flutter(D66, t_mm, b, cr, ct, sweep_deg,
                    h, V_rocket=None, M_rocket=None,
                    do_sweep=True, x_ea=0.40):
    t              = t_mm * 1e-3
    rho, p, T, a   = isa(h)
    if V_rocket is None: V_rocket = M_rocket * a
    else:                M_rocket = V_rocket / a

    AR   = fin_ar(b, cr, ct)
    mac  = fin_mac(b, cr, ct)
    tc   = t / mac
    Geff = 12.0 * D66 / t**3
    sw   = sweep_deg if do_sweep else 0.0

    qr   = flutter_supersonic_quadratic(Geff, rho, t, b, AR, tc, a,
                                        x_ea=x_ea, sweep_deg=sw)
    FSF  = qr['Vf_swept'] / V_rocket

    return dict(h=h, rho=rho, T=T, p=p, a=a,
                V_rocket=V_rocket, M_rocket=M_rocket, q=0.5*rho*V_rocket**2,
                AR=AR, mac=mac, tc=tc, Geff=Geff,
                **qr, FSF=FSF)

# ─── REPORT ─────────────────────────────────────────────────────────────────
L72 = "=" * 72

def report(r, D66, t_mm, cr, ct, b, sweep, fsf_req, do_sweep, x_ea, src):
    print(L72)
    print("  FIN FLUTTER ANALYSIS  [flutterEstimate_v4.py]")
    print(L72)
    if src:
        print(f"  Laminate      : {src}")
    print(f"  D66  = {D66:.4f} N.m     G_eff = {r['Geff']/1e9:.4f} GPa")
    print(f"  t    = {t_mm:.4f} mm")
    print()
    print(f"  FIN GEOMETRY")
    print(f"    cr = {cr*1e3:.1f} mm   ct = {ct*1e3:.1f} mm   b = {b*1e3:.1f} mm")
    print(f"    Λ  = {sweep:.1f} deg   AR = {r['AR']:.4f}   MAC = {r['mac']*1e3:.2f} mm")
    print(f"    t/c = {r['tc']:.4f}  ({r['tc']*100:.2f}%)")
    print(f"    x_EA = {x_ea:.2f}c   e_sub = {x_ea-0.25:.2f}c   e_sup = {0.50-x_ea:.2f}c")
    print()
    print(f"  ISA ATMOSPHERE  h = {r['h']:.1f} m")
    print(f"    ρ = {r['rho']:.4f} kg/m³   T = {r['T']:.2f} K   a = {r['a']:.2f} m/s")
    print(f"    p = {r['p']/1e3:.3f} kPa    q = {r['q']/1e3:.3f} kPa")
    print()
    print(f"  ROCKET  Vr = {r['V_rocket']:.1f} m/s   M = {r['M_rocket']:.4f}")
    print()
    print(f"  SUPERSONIC QUADRATIC  [Ackeret + BAH §5.5]")
    print(f"    S = {r['S']:.1f} m/s   R = {r['R']:.4f}   A = {r['A']:.4f}  "
          f"{'[A > 2 ✓]' if r['A'] > 2 else '[A ≤ 2 ✗]'}")
    print(f"    Δ = {r['disc']:.4f}   u⁺ = {r['u_plus']:.4f}")
    print()
    print(f"    [1] Pre-sweep   : Vf = {r['Vf_super']:.1f} m/s   Mf = {r['Mf_super']:.3f}")
    if do_sweep:
        print(f"    [2] Sweep ×{r['f_sw']:.4f}: Vf = {r['Vf_swept']:.1f} m/s   Mf = {r['Mf_swept']:.3f}")
    else:
        print(f"    [2] Sweep correction : DISABLED")
    print()
    margin = (r['FSF'] / fsf_req - 1.0) * 100
    status = "PASS ✓" if r['FSF'] >= fsf_req else "FAIL ✗"
    print(f"  FLUTTER SAFETY FACTOR")
    print(f"    V_flutter = {r['Vf_swept']:.1f} m/s   V_rocket = {r['V_rocket']:.1f} m/s")
    print(f"    FSF = {r['FSF']:.4f}   [{status}]   required ≥ {fsf_req:.2f}   margin = {margin:+.1f}%")
    print()
    if r['FSF'] >= fsf_req:
        print(f"  ✓  Flutter boundary ({r['Vf_swept']:.0f} m/s) clears flight speed "
              f"({r['V_rocket']:.0f} m/s) by {margin:.1f}%.")
    else:
        print(f"  ✗  Insufficient margin. Target V_flutter ≥ {r['V_rocket']*fsf_req:.0f} m/s")

def altitude_table(D66, t_mm, b, cr, ct, sweep, M0, do_sweep, fsf_req, x_ea):
    print(f"\n{L72}")
    print(f"  ALTITUDE SWEEP  (M_rocket = {M0:.3f} constant)")
    print(L72)
    print(f"  {'h[m]':>6}  {'rho':>8}  {'a[m/s]':>7}  {'Vr[m/s]':>8}  "
          f"{'S[m/s]':>7}  {'A':>6}  {'Vf[m/s]':>8}  {'Mf':>6}  {'FSF':>6}  Status")
    print(f"  {'-'*82}")
    for h in range(0, 4001, 250):
        rho, _, _, a = isa(h)
        Vr = M0 * a
        try:
            r  = compute_flutter(D66, t_mm, b, cr, ct, sweep, h,
                                 V_rocket=Vr, do_sweep=do_sweep, x_ea=x_ea)
            ok = "PASS" if r['FSF'] >= fsf_req else "FAIL"
            tag = " ← max-q" if abs(h - 1462) < 130 else ""
            print(f"  {h:>6d}  {rho:>8.4f}  {a:>7.2f}  {Vr:>8.1f}  "
                  f"{r['S']:>7.1f}  {r['A']:>6.3f}  {r['Vf_swept']:>8.1f}  "
                  f"{r['Mf_swept']:>6.3f}  {r['FSF']:>6.3f}  {ok}{tag}")
        except ValueError:
            print(f"  {h:>6d}  {rho:>8.4f}  {a:>7.2f}  {Vr:>8.1f}  "
                  f"{'—':>7}  {'—':>6}  {'—':>8}  {'—':>6}  {'—':>6}  NO ROOT")

# ─── MAIN ───────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    D66, t_mm, src = None, None, None

    if args.json:
        with open(args.json) as f:
            lam = json.load(f)
        D66  = lam['tailored_beta']['D66_Nm']
        t_mm = lam['tailored_beta']['t_total_mm']
        src  = args.json

    if args.d66 is not None: D66  = args.d66
    if args.t   is not None: t_mm = args.t
    if D66  is None: D66  = 205.96;  print(f"  Default D66 = {D66} N.m")
    if t_mm is None: t_mm = 5.356;   print(f"  Default t   = {t_mm} mm")

    if args.csv:
        max_vel, alt_ft = read_flight_data(args.csv)
        if args.velocity  is None:        args.velocity  = max_vel
        if args.altitude  == 1462.0:      args.altitude  = alt_ft * 0.3048

    do_sweep = not args.sweep_off
    sweep    = 0.0 if args.sweep_off else args.sweep

    r = compute_flutter(D66, t_mm, args.span, args.cr, args.ct, sweep,
                        args.altitude,
                        V_rocket=args.velocity, M_rocket=args.mach,
                        do_sweep=do_sweep, x_ea=args.x_ea)
    print()
    report(r, D66, t_mm, args.cr, args.ct, args.span, sweep,
           args.fsf_req, do_sweep, args.x_ea, src)

    if args.sweep_table:
        altitude_table(D66, t_mm, args.span, args.cr, args.ct, sweep,
                       r['M_rocket'], do_sweep, args.fsf_req, args.x_ea)

    print(f"\n{L72}")

if __name__ == "__main__":
    main()
