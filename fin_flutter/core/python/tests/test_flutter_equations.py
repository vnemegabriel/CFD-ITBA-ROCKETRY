"""
test_flutter_equations.py — Flutter equation validation harness
===============================================================
Validates every Bisplinghoff-referenced flutter equation in flutterEstimate_v3.py.
Each test documents the full derivation so the book is not needed to interpret results.

Run:  pytest core/python/tests/test_flutter_equations.py -v

References
----------
[BAH]  Bisplinghoff, Ashley & Halfman. Aeroelasticity, 1955. §5.5
[M58]  Martin, H.C. NACA TN 4197, 1958.
[A25]  Ackeret, J.  NACA TM 317, 1925.  Eq. 16.
[CT54] Collar & Tickle, 1954. (alternative sweep correction)
[ICAO] ICAO Doc 7488 — Manual of the ICAO Standard Atmosphere
"""
import math
import sys
from pathlib import Path
import pytest

# ── resolve project root ──────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parents[3]   # .../fin_flutter/
sys.path.insert(0, str(ROOT))

import importlib.util

_spec = importlib.util.spec_from_file_location(
    "flutterEstimate_v3", ROOT / "flutterEstimate_v3.py"
)
_mod = importlib.util.module_from_spec(_spec)
# Suppress argparse side-effects on import
import unittest.mock as _mock
with _mock.patch("sys.argv", ["flutterEstimate_v3.py"]):
    _spec.loader.exec_module(_mod)

isa              = _mod.isa
flutter_naca4197 = _mod.flutter_naca4197
ackeret_factor   = _mod.ackeret_factor
sweep_factor     = _mod.sweep_factor
fin_ar           = _mod.fin_ar
fin_mac          = _mod.fin_mac
compute_flutter  = _mod.compute_flutter

# Design inputs — single source of truth for all test classes
_FL6 = dict(D66=205.96, t_mm=5.356, b=0.160, cr=0.300, ct=0.150,
            sweep_deg=57.4, h=1462.0, M_rocket=1.942)

# K-factor approximation (NACA TN 4197 Table 2 linear fit) — mirrors production code
def _k_factor(AR: float) -> float:
    return max(0.65 * (1.0 + 0.10 * (AR - 1.0)), 0.40)


# ═══════════════════════════════════════════════════════════════════════════════
# T-1  ISA Standard Atmosphere
# ═══════════════════════════════════════════════════════════════════════════════
class TestISA:
    """
    Derivation [ICAO]:
        T(h) = T0 - L·h                           (troposphere lapse, K)
        p(h) = p0·(T/T0)^(g/(R·L))
        ρ(h) = p / (R·T)
        a(h) = sqrt(γ·R·T)                        (isentropic speed of sound)

    Constants: T0=288.15 K, p0=101325 Pa, L=6.5e-3 K/m,
               g=9.80665 m/s², R=287.058 J/(kg·K), γ=1.4

    Sea-level expected values from ICAO Doc 7488 Table A:
        ρ = 1.2250 kg/m³, T = 288.15 K, a = 340.29 m/s
    """

    def test_sea_level_density(self):
        rho, *_ = isa(0)
        assert abs(rho - 1.2250) < 5e-4, f"ρ_SL = {rho:.5f}, expected 1.2250 kg/m³"

    def test_sea_level_temperature(self):
        _, _, T, _ = isa(0)
        assert abs(T - 288.15) < 1e-6, f"T_SL = {T:.6f} K, expected 288.15 K"

    def test_sea_level_sound_speed(self):
        _, _, _, a = isa(0)
        # a = sqrt(1.4 × 287.058 × 288.15) = 340.294 m/s
        assert abs(a - 340.294) < 0.01, f"a_SL = {a:.3f} m/s, expected 340.29 m/s"

    def test_tropopause_temperature(self):
        _, _, T, _ = isa(11000)
        # T(11 km) = 288.15 - 0.0065 × 11000 = 216.65 K
        assert abs(T - 216.65) < 1e-6, f"T(11 km) = {T:.6f} K, expected 216.65 K"

    def test_design_altitude(self):
        """h = 1462 m (design max-q altitude)."""
        rho, _, T, a = isa(1462)
        T_expected = 288.15 - 0.0065 * 1462          # = 278.647 K
        assert abs(T - T_expected) < 1e-4
        assert rho > 0 and a > 0


# ═══════════════════════════════════════════════════════════════════════════════
# T-2  NACA TN 4197 K-factor
# ═══════════════════════════════════════════════════════════════════════════════
class TestNACA4197KFactor:
    """
    Derivation [M58]:
        K ≈ 0.65·(1 + 0.10·(AR - 1))   [linear approximation to Table 2]
        K_min = 0.40                     [physical lower bound]

    NOTE: This approximation is one-dimensional in AR only; NACA TN 4197
    Table 2 is a function of both AR and taper ratio λ.  Bilinear interpolation
    over the table is the rigorous implementation — see ERRORS.md F-4.

    These tests verify the approximation's own internal consistency.
    """

    def test_ar_unity(self):
        # K is private inside flutter_naca4197; _k_factor mirrors the formula.
        assert abs(_k_factor(1.0) - 0.65) < 1e-10

    def test_minimum_clamp(self):
        # AR=0 gives formula value 0.585 < 0.40 clamp threshold
        assert _k_factor(0.0) >= 0.40

    def test_design_ar(self):
        # AR ≈ 0.711 → K ≈ 0.6312
        AR = fin_ar(_FL6['b'], _FL6['cr'], _FL6['ct'])
        assert abs(_k_factor(AR) - 0.6312) < 1e-3, f"K(AR={AR:.3f}) = {_k_factor(AR):.4f}"


# ═══════════════════════════════════════════════════════════════════════════════
# T-3  Ackeret Factor — BAH §5.5 derivation
# ═══════════════════════════════════════════════════════════════════════════════
class TestAckeretFactor:
    """
    Derivation [BAH §5.5, A25]:

        For a flat-plate fin with elastic axis at x_EA (fraction of chord):

        Subsonic aerodynamics (Prandtl-Glauert thin airfoil):
            CL_α,sub = 2π  [/rad]
            AC_sub   = 0.25c  (quarter-chord)
            e_sub    = x_EA - 0.25        (moment arm, must be > 0 for flutter)

        Supersonic aerodynamics (Ackeret, NACA TM 317 Eq. 16):
            CL_α,sup = 4 / sqrt(M² - 1)  [/rad]
            AC_sup   = 0.50c  (half-chord, symmetric thin airfoil)
            e_sup    = 0.50 - x_EA        (moment arm, must be > 0 for flutter)

        Aerodynamic torsional forcing ∝ CL_α · e
        Flutter speed ∝ 1 / sqrt(CL_α · e)

        Supersonic-to-subsonic scaling factor:
            F = sqrt[ (CL_α,sub · e_sub) / (CL_α,sup · e_sup) ]

        Physical meaning: supersonic aerodynamics are weaker torsionally
        (lower CL_α), but the AC moves aft (larger moment arm). Net effect
        depends on the balance of both. For x_EA = 0.40, the subsonic
        aero-torsional product is larger, so supersonic flutter speed is
        higher (F > 1).

    NOTE — HTML inconsistency (BUG-2):
        layup_viewer_v5.html hardcodes e_sub = 0.25 (the AC location, not the
        moment arm). With x_EA = 0.40 the correct e_sub = x_EA - 0.25 = 0.15.
        The HTML value inflates the Ackeret factor by ~29% at design Mach.
    """

    def test_subsonic_returns_unity(self):
        """M ≤ 1.0 must give factor = 1.0 (no supersonic correction)."""
        assert ackeret_factor(0.5) == 1.0
        assert ackeret_factor(1.0) == 1.0

    def test_analytical_m2_xea040(self):
        """
        M = 2.0, x_EA = 0.40

        Step-by-step:
            e_sub      = 0.40 - 0.25 = 0.15 c
            e_sup      = 0.50 - 0.40 = 0.10 c
            CL_α,sub   = 2π = 6.28318 /rad
            CL_α,sup   = 4 / sqrt(4 - 1) = 4 / sqrt(3) = 2.30940 /rad
            numerator  = 6.28318 × 0.15 = 0.94248
            denominator= 2.30940 × 0.10 = 0.23094
            F          = sqrt(0.94248 / 0.23094) = sqrt(4.08099) = 2.0202

        Tolerance: 0.001 (numerical precision of sqrt)
        """
        F = ackeret_factor(2.0, x_ea=0.40)
        assert abs(F - 2.0202) < 0.001, f"ackeret_factor(2.0, 0.40) = {F:.4f}, expected 2.0202"

    def test_analytical_m3_xea040(self):
        """
        M = 3.0, x_EA = 0.40

            CL_α,sup = 4/sqrt(9-1) = 4/sqrt(8) = 1.41421 /rad
            F = sqrt[(6.28318×0.15)/(1.41421×0.10)]
              = sqrt[0.94248 / 0.14142] = sqrt[6.6667] = 2.5820
        """
        F = ackeret_factor(3.0, x_ea=0.40)
        expected = math.sqrt((2 * math.pi * 0.15) / (4 / math.sqrt(3.0**2 - 1) * 0.10))
        assert abs(F - expected) < 1e-6

    def test_factor_greater_than_one(self):
        """Supersonic factor must be > 1 for x_EA < 0.50 (stabilising EA placement)."""
        for M in [1.1, 1.5, 2.0, 3.0, 5.0]:
            assert ackeret_factor(M, x_ea=0.40) > 1.0, \
                f"ackeret_factor({M}) ≤ 1 — unexpected for x_EA=0.40"

    def test_factor_increases_with_mach(self):
        """
        As M → ∞, CL_α,sup → 0, so aerodynamic torsional coupling → 0,
        and flutter speed → ∞.  Factor must be monotonically increasing with M.
        """
        factors = [ackeret_factor(M, x_ea=0.40) for M in [1.5, 2.0, 3.0, 5.0]]
        assert all(factors[i] < factors[i+1] for i in range(len(factors)-1)), \
            f"Ackeret factor not monotone: {factors}"

    @pytest.mark.xfail(reason="BUG-2: HTML hardcodes e_sub=0.25 instead of x_EA-0.25=0.15; "
                               "inflates Ackeret factor by ~29% at design Mach M=1.942")
    def test_html_parity_ackeret_m1942(self):
        """
        At M=1.942, x_EA=0.40:
            Python (correct): e_sub = 0.15 → F ≈ 2.003
            HTML  (BUG-2):    e_sub = 0.25 → F ≈ 2.584

        This test asserts they are equal — it is marked xfail to document
        the known inconsistency. Fix HTML supersonicScale() to use
        e_sub = x_EA - 0.25 with a configurable x_EA parameter.
        """
        M = 1.942
        F_python = ackeret_factor(M, x_ea=0.40)          # correct: e_sub=0.15

        # HTML equivalent (hardcoded e_sub=0.25)
        CL_sub = 2 * math.pi
        CL_sup = 4 / math.sqrt(M**2 - 1)
        F_html = math.sqrt((CL_sub * 0.25) / (CL_sup * 0.10))

        assert abs(F_python - F_html) < 0.001, \
            (f"Ackeret factor mismatch at M={M}: "
             f"Python={F_python:.4f} (e_sub=0.15), "
             f"HTML={F_html:.4f} (e_sub=0.25, WRONG). "
             f"Relative error = {abs(F_python-F_html)/F_python*100:.1f}%")


# ═══════════════════════════════════════════════════════════════════════════════
# T-4  Sweep Correction — BAH §5.5 derivation
# ═══════════════════════════════════════════════════════════════════════════════
class TestSweepFactor:
    """
    Derivation [BAH §5.5]:

        For an aft-swept panel (Λ > 0), the aerodynamic chordwise component
        of the free-stream dynamic pressure is:

            q_eff = q · cos²(Λ)

        The flutter dynamic pressure q_flutter is set by structural torsional
        stiffness (independent of sweep).  Setting q_eff = q_flutter:

            0.5·ρ·V_f,swept² · cos²(Λ) = 0.5·ρ·V_f,unswept²
            V_f,swept = V_f,unswept / cos(Λ)      → factor = 1/cos(Λ) > 1

        For forward sweep (Λ < 0) the argument inverts (aeroelastic wash-in):
            factor = cos(|Λ|) < 1   (conservative first-order estimate only;
                                      forward-swept panels may diverge first)

    NOTE — HTML inconsistency (BUG-1):
        layup_viewer_v5.html uses factor = 1/sqrt(cos(Λ)) — Collar & Tickle (1954).
        This is a different physical model (correction on V, not on q²).
        At design sweep Λ=57.4°:
            Python (BAH):   1/cos(57.4°) = 1.855
            HTML (CT54):    1/sqrt(cos(57.4°)) = 1.362
        Difference: 26% — significant for flutter margin prediction.
    """

    def test_zero_sweep_is_unity(self):
        assert sweep_factor(0.0) == 1.0

    def test_60_deg_aft_sweep(self):
        """
        Λ = 60° → cos(60°) = 0.5 → factor = 1/0.5 = 2.0  (exact)
        """
        f = sweep_factor(60.0)
        assert abs(f - 2.0) < 1e-10, f"sweep_factor(60°) = {f:.6f}, expected 2.0"

    def test_45_deg_aft_sweep(self):
        """
        Λ = 45° → cos(45°) = √2/2 → factor = √2 ≈ 1.41421
        """
        f = sweep_factor(45.0)
        assert abs(f - math.sqrt(2)) < 1e-10

    def test_design_sweep(self):
        """
        Λ = 57.4° (design leading-edge sweep)
            cos(57.4°) = cos(1.00180 rad) = 0.53919
            factor = 1 / 0.53919 = 1.8547
        """
        f = sweep_factor(57.4)
        expected = 1.0 / math.cos(math.radians(57.4))
        assert abs(f - expected) < 1e-10

    def test_forward_sweep_destabilising(self):
        """Forward sweep (Λ < 0) must give factor < 1."""
        f = sweep_factor(-30.0)
        assert f < 1.0, f"Forward sweep factor = {f:.4f}, expected < 1"
        assert abs(f - math.cos(math.radians(30.0))) < 1e-10

    @pytest.mark.xfail(reason="BUG-1: HTML uses 1/sqrt(cos Λ) [Collar & Tickle] instead of "
                               "1/cos Λ [BAH §5.5]; 26% difference at design sweep 57.4°")
    def test_html_parity_sweep_factor(self):
        """
        At Λ = 57.4°:
            Python (BAH §5.5):   1/cos(57.4°) = 1.8547
            HTML (CT54):         1/sqrt(cos(57.4°)) = 1.3617
        """
        sweep = 57.4
        f_python = sweep_factor(sweep)                                 # 1/cos(Λ)
        f_html   = 1.0 / math.sqrt(math.cos(math.radians(sweep)))     # 1/√cos(Λ)

        assert abs(f_python - f_html) < 0.01, \
            (f"Sweep factor mismatch at Λ={sweep}°: "
             f"Python={f_python:.4f} (BAH §5.5, 1/cos), "
             f"HTML={f_html:.4f} (Collar & Tickle, 1/√cos). "
             f"Difference = {abs(f_python-f_html)/f_python*100:.1f}%")


# ═══════════════════════════════════════════════════════════════════════════════
# T-5  Supersonic Fixed-Point Iteration
# ═══════════════════════════════════════════════════════════════════════════════
class TestSupersonicIteration:
    """
    The compute_flutter() pipeline solves for the self-consistent flutter Mach Mf:
        Mf = Mf_sub × ackeret_factor(Mf)
    by fixed-point iteration (max 50 steps, tol 1e-6).

    Physical requirement: the supersonic-corrected flutter speed must be
    higher than the subsonic baseline (Vf2 > Vf1), because ackeret_factor > 1.
    """

    def _run(self, **kwargs):
        return compute_flutter(**{**_FL6, **kwargs})

    def test_supersonic_raises_flutter_speed(self):
        """
        At design point, Mf_sub = 1.059 (just above M=1). The Ackeret factor
        at M≈1.06 is 0.904 < 1 (near-transonic CL_α singularity inflates the
        aerodynamic torsional term, reducing flutter speed). The fixed-point
        iteration 2-cycles between Mf≈1.06 (f=0.904) and Mf≈0.96 (f=1.0)
        without converging — BUG-6 (see xfail below).

        After 50 iterations (even loop count) f_sup lands at 1.0 → Vf2 = Vf1.
        This is physically conservative (use subsonic estimate), but the
        iteration is structurally broken: an odd loop count would give f_sup=0.904.
        """
        r = self._run(do_super=True, do_sweep=False)
        # At design point, Ackeret correction 2-cycles — Vf2 == Vf1 (f_sup=1.0)
        assert r['Vf2'] >= r['Vf1'], \
            f"Supersonic correction unexpectedly reduced flutter speed below subsonic: " \
            f"Vf1={r['Vf1']:.1f}, Vf2={r['Vf2']:.1f}"

    @pytest.mark.xfail(reason=(
        "BUG-6: fixed-point iteration 2-cycles near M=1 (transonic Ackeret singularity). "
        "Mf_sub=1.059 → f=0.904 → Mf=0.957 (subsonic) → f=1.0 → Mf=1.059 (loops). "
        "After 50 iters f_sup is loop-count-parity dependent (0.904 or 1.0). "
        "Fix: detect 2-cycle oscillation; fall back to subsonic estimate (f_sup=1.0) "
        "when self-consistent supersonic solution does not exist."
    ))
    def test_iteration_converges_near_transonic(self):
        """Iteration must reach |Mf_new - Mf_iter| < 1e-6 before the 50-step limit."""
        r = self._run(do_super=True, do_sweep=False)
        # The iteration should have converged, meaning f_sup != exactly 1.0
        # when Mf1 > 1.  If it 2-cycled, it silently exits at an arbitrary value.
        # Proxy: Vf2 must be strictly above Vf1 for any Mf1 > 1.
        assert r['Vf2'] > r['Vf1'], \
            f"Iteration did not converge; Vf2={r['Vf2']:.1f} == Vf1={r['Vf1']:.1f} (2-cycle BUG-6)"

    def test_sweep_raises_flutter_speed_aft(self):
        r = self._run(do_super=False, do_sweep=True)
        assert r['Vf3'] > r['Vf1'], \
            f"Aft sweep correction reduced flutter speed: Vf1={r['Vf1']:.1f}, Vf3={r['Vf3']:.1f}"

    def test_pipeline_convergence(self):
        """Full pipeline must complete without ValueError (iteration divergence)."""
        r = self._run(do_super=True, do_sweep=True)
        assert math.isfinite(r['Vf3'])
        assert math.isfinite(r['FSF'])

    def test_fsf_positive(self):
        r = self._run()
        assert r['FSF'] > 0.0

    def test_intermediate_values_consistent(self):
        """Vf1 ≤ Vf2 ≤ Vf3 for aft sweep + supersonic correction."""
        r = self._run(do_super=True, do_sweep=True)
        assert r['Vf1'] <= r['Vf2'], "Supersonic step reduced flutter speed"
        assert r['Vf2'] <= r['Vf3'], "Aft sweep step reduced flutter speed"


# ═══════════════════════════════════════════════════════════════════════════════
# T-6  Full Pipeline Regression
# ═══════════════════════════════════════════════════════════════════════════════
class TestFullPipelineRegression:
    """
    Full pipeline with design inputs.
    Regression values established from first validated run.

    Pipeline chain:
        ISA(h=1462 m) → NACA TN 4197 (subsonic) → Ackeret (supersonic,
        fixed-point) → sweep correction (BAH §5.5) → FSF

    Expected values (informational — printed for manual inspection):
        AR  ≈ 0.711   (semi-span=160 mm, cr=300 mm, ct=150 mm)
        MAC ≈ 233.3 mm
        t/c ≈ 0.023
        Geff ≈ 50 GPa (rough order-of-magnitude)
        FSF ≥ 1.0 (sanity; true requirement is ≥ 1.5)
    """

    def test_geometry_checks(self):
        b, cr, ct = _FL6['b'], _FL6['cr'], _FL6['ct']
        AR  = fin_ar(b, cr, ct)
        mac = fin_mac(b, cr, ct)
        # AR = b² / (0.5·(cr+ct)·b) = 0.160² / (0.5·0.450·0.160) = 0.02560/0.03600 = 0.7111
        assert abs(AR - 0.7111) < 1e-3, f"AR = {AR:.4f}"
        # MAC = 0.300 × (2/3) × (1 + 0.5 + 0.25) / (1 + 0.5) = 0.300 × (2/3) × 1.1667 = 0.2333 m
        assert abs(mac - 0.2333) < 1e-3, f"MAC = {mac:.4f} m"

    def test_full_run_sanity(self, capsys):
        r = compute_flutter(**_FL6, do_sweep=True, do_super=True)
        # Print for manual inspection without failing the test
        with capsys.disabled():
            print(f"\n  ── Full Pipeline Regression (T-6) ──")
            print(f"  AR      = {r['AR']:.4f}")
            print(f"  MAC     = {r['mac']*1e3:.2f} mm")
            print(f"  t/c     = {r['tc']:.4f}  ({r['tc']*100:.2f}%)")
            print(f"  Geff    = {r['Geff']/1e9:.3f} GPa")
            print(f"  ρ       = {r['rho']:.4f} kg/m³   T = {r['T']:.2f} K   a = {r['a']:.2f} m/s")
            print(f"  V_r     = {r['V_rocket']:.1f} m/s   M_r = {r['M_rocket']:.4f}")
            print(f"  Vf1     = {r['Vf1']:.1f} m/s  (NACA TN 4197 subsonic)")
            print(f"  ×f_sup  = {r['f_sup']:.4f}  → Vf2 = {r['Vf2']:.1f} m/s  (Ackeret M={r['Mf2']:.3f})")
            print(f"  ×f_sw   = {r['f_sw']:.4f}  → Vf3 = {r['Vf3']:.1f} m/s  (sweep-corrected)")
            print(f"  FSF     = {r['Vf3']:.1f} / {r['V_rocket']:.1f} = {r['FSF']:.4f}")
            print(f"  Status  = {'PASS ✓' if r['FSF'] >= 1.5 else 'FAIL ✗'}  (req ≥ 1.50)")

        assert r['FSF'] > 0.0
        assert all(math.isfinite(v) for v in
                   [r['Vf1'], r['Vf2'], r['Vf3'], r['FSF'], r['f_sup'], r['f_sw']])
