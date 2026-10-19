"""Energy per cycle at finite amplitude: aerodynamic work against structural dissipation, and the LCO amplitude."""
import re
from pathlib import Path

import numpy as np

from . import cfd


def signal(run):
    """amplitude, omega, coef of the harmonic signal in <run>/constant/flutterModes."""
    text = (Path(run) / "constant" / "flutterModes").read_text()
    sig = text[text.index("signal"):]
    get = lambda k: re.search(rf"^\s*{k}\s+(.+?);", sig, re.M).group(1)
    if get("type") != "harmonic":
        raise ValueError(f"{run}: the signal is not harmonic")
    coef = np.array(get("coef").strip("()").split(), float)
    omegas = np.array(re.search(r"^omega\s+\((.+?)\);", text, re.M).group(1).split(), float)
    return float(get("amplitude")), float(get("omega")), coef, omegas


def cycle_work(run, n_cycles=2):
    """Work the air does on the fin per cycle, W_a = (1/N) int sum_k F_k deta_k/dt dt over the last N cycles."""
    A, w, coef, omegas = signal(run)
    t, eta, F = cfd.read_genforces(run)
    T = 2 * np.pi / w
    if t[-1] - t[0] < (n_cycles + 1) * T:
        raise ValueError(f"{run}: shorter than {n_cycles + 1} cycles")
    sel = t >= t[-1] - n_cycles * T
    etad = np.gradient(eta[sel], t[sel], axis=0)
    W = np.trapezoid(np.sum(F[sel] * etad, axis=1), t[sel]) / n_cycles
    return {"run": str(run), "A": A, "omega": w, "coef": coef, "omegas": omegas, "W_a": W}


def dissipation(coef, omegas, g):
    """Energy per cycle of a hysteretic loss factor g_k: pi sum g_k Omega_k^2 |eta_k|^2 (mass-normalised modes)."""
    return np.pi * np.sum(np.asarray(g) * omegas[: len(coef)] ** 2 * coef ** 2)


def balance(runs, g, q_scale=1.0):
    """Table over amplitude: W_a (times q_scale = q_flight/q_cfd), W_d, ratio; LCO amplitudes where they cross."""
    rows = sorted((cycle_work(r) for r in runs), key=lambda d: d["A"])
    for d in rows:
        d["W_a"] *= q_scale
        d["W_d"] = dissipation(d["coef"], d["omegas"], g)
        d["ratio"] = d["W_a"] / d["W_d"]
    lco = []
    for a, b in zip(rows, rows[1:]):
        fa, fb = a["W_a"] - a["W_d"], b["W_a"] - b["W_d"]
        if fa * fb < 0:
            A = a["A"] + (b["A"] - a["A"]) * fa / (fa - fb)
            lco.append({"A": A, "stable": fb < fa})
    return rows, lco
