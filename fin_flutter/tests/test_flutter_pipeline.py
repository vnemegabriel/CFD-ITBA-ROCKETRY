import re
from pathlib import Path

import numpy as np
import pytest
from scipy.io import wavfile

from finflutter import aero, aeroelastic, cfd, gaf, geometry, lco, static, taptest
from finflutter.run import build_structure, load_case

ROOT = Path(__file__).parents[1]
MESHINFO = ROOT.parent / "flightModel/mesh/meshInfo"


@pytest.fixture(scope="module")
def s():
    c = load_case(ROOT / "cases/aconcagua.yaml")
    c["mesh"], c["modes"] = {"nx": 20, "ny": 8}, 4
    return build_structure(c)


def test_fin_matches_flightmodel():
    fin = geometry.fin_from_flightmodel()
    x0 = float(re.search(r"^finX0\s+([0-9.]+);", MESHINFO.read_text(), re.M).group(1))
    assert fin.x_root_le == pytest.approx(x0, abs=1e-5)
    assert fin.sweep_deg == pytest.approx(57.46, abs=0.01)
    for eta in (0, 1):
        c = fin.chord(eta)
        assert 2 * fin.half_thickness(c / 2, c) == pytest.approx(fin.t)
        assert 2 * fin.half_thickness(fin.rho, c) == pytest.approx(2 * fin.rho, rel=0.01)


def test_mode_interpolation_reproduces_fe_nodes(s):
    fin, nodes = s["fin"], s["nodes"]
    r = np.radians(45.0)
    er, en = np.array([0, np.cos(r), np.sin(r)]), np.array([0, -np.sin(r), np.cos(r)])
    pts = (fin.x_root_le + nodes[:, 0])[:, None] * np.array([1, 0, 0]) + (fin.R + nodes[:, 1])[:, None] * er + 0.003 * en
    w, _ = static.shape(pts, s, 45.0)
    assert np.allclose(w, s["phi"], atol=1e-9 * np.abs(s["phi"]).max())


def _fake_run(tmp, t, eta, F, M=1.5, rho=1.0, U=500.0, signal=None):
    run = tmp / f"run{len(list(tmp.iterdir()))}"
    (run / "postProcessing/genForces/0").mkdir(parents=True)
    (run / "constant").mkdir()
    (run / "run.info").write_text(f"case x;\nflow\n{{\n    Minf {M};\n    rhoInf {rho};\n    Uinf {U};\n    pInf 1e5;\n    alpha 0;\n}}\n")
    np.savetxt(run / "postProcessing/genForces/0/genForces.dat", np.column_stack([t, eta, F]), header="t eta F")
    if signal:
        (run / "constant/flutterModes").write_text(signal)
    return run


def test_gaf_identification_recovers_known_q(tmp_path):
    rng = np.random.default_rng(1)
    n, b, U, rho = 2, 0.1, 500.0, 1.0
    q = 0.5 * rho * U ** 2
    A0, A1 = rng.normal(size=(n, n)), rng.normal(size=(n, n))
    sig = cfd.pulse(n, 1, 5e-4, np.ones(n), 1500.0)
    t = np.arange(0, 0.02, 1e-6)
    runs = []
    for j in range(n):
        eta = np.zeros((len(t), n))
        eta[:, j] = 1e-4 * np.exp(-((t - sig["t0"]) / sig["tau"]) ** 2)
        etad = np.gradient(eta, t, axis=0)
        F = -q * (eta @ A0.T + (b / U) * etad @ A1.T)
        runs.append(_fake_run(tmp_path, t, eta, F, rho=rho, U=U))
    d = gaf.assemble(runs, b=b)
    exact = A0[None] + 1j * d["k"][:, None, None] * A1[None]
    assert np.abs(d["Q"] - exact).max() < 1e-3 * np.abs(exact).max()
    fit = aeroelastic.roger_fit(d["k"], d["Q"], [0.1, 0.3])
    assert np.allclose(fit["A0"], A0, atol=1e-3) and np.allclose(fit["A1"], A1, atol=1e-3)


def test_roger_state_space_equals_exact_piston_roots(s):
    om, g = s["omega"], 0.02 * np.ones(len(s["omega"]))
    pt, b = aero.PistonTheory(s["phi"], s["mats"]), 0.1125
    gafs = [{"mach": M, "b": b, "k": np.linspace(0, 1, 20),
             "Q": np.array([pt(M)[0] + 1j * k * pt(M)[1] / b for k in np.linspace(0, 1, 20)])} for M in (1.6, 2.0)]
    rfa = aeroelastic.RogerTable(gafs)(2.0)
    r1 = np.sort_complex(aeroelastic.roots(om, g, *pt(2.0), 2e5, 600.0))
    r2 = np.linalg.eigvals(aeroelastic.state_matrix(om, g, rfa, 2e5, 600.0))
    r2 = np.sort_complex(r2[np.abs(r2.imag) > 1.0])
    assert np.allclose(r1, r2, rtol=1e-8)


def test_tap_test_recovers_frequency_and_loss_factor(tmp_path):
    fs = 44100
    t = np.arange(0, 0.5, 1 / fs)
    x = np.zeros(int(2.5 * fs))
    for t0 in (0.2, 1.2):
        sgl = sum(a * np.exp(-z * 2 * np.pi * f * t) * np.sin(2 * np.pi * f * t)
                  for a, f, z in ((1, 240, 0.01), (0.6, 800, 0.015)))
        x[int(t0 * fs): int(t0 * fs) + len(t)] += sgl
    wavfile.write(tmp_path / "tap.wav", fs, (x / np.abs(x).max() * 30000).astype(np.int16))
    m = taptest.summary(taptest.analyse(tmp_path / "tap.wav", n_modes=2), [237, 795])
    assert m[0]["f_mean"] == pytest.approx(240, abs=2) and m[0]["g_mean"] == pytest.approx(0.02, abs=0.002)
    assert m[1]["f_mean"] == pytest.approx(800, abs=4) and m[1]["g_mean"] == pytest.approx(0.03, abs=0.003)


def test_energy_per_cycle_of_a_damping_force(tmp_path):
    w, c, b, U, q, a1 = 2 * np.pi * 300, 1e-4, 0.1, 500.0, 1e5, 0.4
    t = np.arange(0, 6 * 2 * np.pi / w, 1e-6)
    eta = (c * np.sin(w * t))[:, None]
    F = -q * a1 * (b / U) * (c * w * np.cos(w * t))[:, None]
    sig = f"omega ({w});\nsignal\n{{\n    type harmonic;\n    amplitude 5e-4;\n    coef ({c});\n    omega {w};\n}}\n"
    run = _fake_run(tmp_path, t, eta, F, signal=sig)
    W = lco.cycle_work(run)["W_a"]
    assert W == pytest.approx(-q * a1 * (b / U) * c ** 2 * w * np.pi, rel=1e-3)
