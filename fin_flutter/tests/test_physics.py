from pathlib import Path

import numpy as np
import pytest

from finflutter import aero, aeroelastic, laminate, structure
from finflutter.run import run

ROOT = Path(__file__).parents[1]
ISO = laminate.isotropic(70e9, 0.3, 0.005, 2700)
D = ISO["D"][0, 0]


def square(n):
    nodes, elems = structure.mesh(1, 1, 1, 0, n, n)
    return nodes, structure.assemble(nodes, elems, ISO)


def test_leissa_cantilever_plate():
    nodes, m = square(16)
    w, _ = structure.modes(m, structure.clamped_dofs(nodes, [[0, 1]]), 5)
    lam = w * np.sqrt(ISO["rho"] * ISO["t"] / D)
    assert np.allclose(lam, [3.471, 8.507, 21.29, 27.20, 30.96], rtol=0.015)


def test_dowell_simply_supported_panel_flutter():
    nodes, m = square(16)
    x, y = nodes.T
    edge = np.flatnonzero((x < 1e-9) | (x > 1 - 1e-9) | (y < 1e-9) | (y > 1 - 1e-9))
    w, phi = structure.modes(m, 3 * edge, 12)
    A0 = 2 * phi.T @ (m["S"] @ phi)  # one-sided, beta = 1
    q, kind, _ = aeroelastic.critical_q(w, 0 * w, A0, 0 * A0, 1.0, 1000 * D)
    assert kind == "flutter" and 2 * q / D == pytest.approx(512, rel=0.03)


def test_piston_theory_damps_single_mode():
    nodes, m = square(8)
    w, phi = structure.modes(m, structure.clamped_dofs(nodes, [[0, 1]]), 1)
    A0, B1 = aero.PistonTheory(phi, m)(2.0)
    re = [aeroelastic.roots(w, 0.02 * w, A0, B1, q, 600.0).real.max() for q in (0, 1e4, 1e5)]
    assert re[0] > re[1] > re[2]


def test_tabulated_gaf_recovers_rational_form():
    nodes, m = square(8)
    _, phi = structure.modes(m, structure.clamped_dofs(nodes, [[0, 1]]), 3)
    pt, b, ks, machs = aero.PistonTheory(phi, m), 0.5, [0.0, 0.1, 0.3], [1.5, 2.0]
    Q = [[pt(M)[0] + 1j * k * pt(M)[1] / b for k in ks] for M in machs]
    A0, B1 = aero.TabulatedGAF(machs, ks, Q, b)(1.75)
    assert np.allclose(A0, (pt(1.5)[0] + pt(2.0)[0]) / 2) and np.allclose(B1, (pt(1.5)[1] + pt(2.0)[1]) / 2)


def test_laminate_matches_reference_clt():
    L = laminate.laminate("ar1", 8.0, 0.5, 5.0)
    assert np.allclose(L["D"][[0, 1, 2, 0, 1], [0, 1, 2, 2, 2]], [2018.418, 2039.921, 656.716, 41.336, 19.638], rtol=1e-4)
    assert L["rho"] == pytest.approx(1501.17, abs=0.01)


def test_pipeline_runs():
    import yaml
    c = yaml.safe_load((ROOT / "cases/ar1_8mm.yaml").read_text())
    c["mesh"] = {"nx": 12, "ny": 6}
    c["modes"] = 4
    r = run(c, ROOT)
    assert r["f_Hz"][0] == pytest.approx(237, rel=0.05) and len(r["envelope"]) > 100
