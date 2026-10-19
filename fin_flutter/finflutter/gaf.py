"""Generalised aerodynamic forces from forced-motion runs: Q_ij(k) = -F_i(w) / (q eta_j(w))."""
from pathlib import Path

import numpy as np

from . import cfd


def excited_mode(eta):
    amp = np.abs(eta).max(axis=0)
    if amp.max() == 0:
        return None
    if np.sum(amp > 1e-9 * amp.max()) > 1:
        raise ValueError("more than one mode moves: the pulse identification needs one mode per run")
    return int(np.argmax(amp))


def response(run, zero=None, dt=None):
    """Uniform-time (t, eta, F) of a run, minus the A = 0 run (or minus F at the first time)."""
    t, eta, F = cfd.read_genforces(run)
    if zero is not None:
        t0, _, F0 = cfd.read_genforces(zero)
        if t0[-1] < t[-1] - 1e-12:
            raise ValueError(f"{zero} ends at t = {t0[-1]:g}, before {run} (t = {t[-1]:g})")
        F = F - np.column_stack([np.interp(t, t0, F0[:, i]) for i in range(F.shape[1])])
    else:
        F = F - F[0]
    dt = dt or float(np.median(np.diff(t)))
    tu = np.arange(t[0], t[-1], dt)
    return tu, np.column_stack([np.interp(tu, t, c) for c in eta.T]), np.column_stack([np.interp(tu, t, c) for c in F.T])


def column(run, zero=None, frac=0.01, pad=4):
    """One column of Q from a pulse run: (j, f [Hz], Q[:, j](f)) where the pulse spectrum is above frac of its peak."""
    t, eta, F = response(run, zero)
    j = excited_mode(eta)
    if j is None:
        raise ValueError(f"{run}: no mode moves (signal none)")
    flow = cfd.run_flow(run)
    q = 0.5 * flow["rhoInf"] * flow["Uinf"] ** 2
    n = pad * len(t)
    dt = t[1] - t[0]
    E = np.fft.rfft(eta[:, j], n)
    Fh = np.fft.rfft(F, n, axis=0)
    f = np.fft.rfftfreq(n, dt)
    keep = np.abs(E) > frac * np.abs(E).max()
    return j, f[keep], -Fh[keep] / (q * E[keep, None]), flow


def assemble(runs, zero=None, b=None, nk=40):
    """Q[k, i, j] on a common k grid from one pulse run per mode. Returns a dict for np.savez."""
    cols, flow = {}, None
    for r in runs:
        j, f, Qj, flow = column(r, zero)
        if j in cols:
            raise ValueError(f"mode {j+1} moves in two runs: give the linearity runs to 'linearity', not here")
        cols[j] = (f, Qj)
    n = next(iter(cols.values()))[1].shape[1]
    missing = sorted(set(range(n)) - set(cols))
    U = flow["Uinf"]
    fmax = min(c[0][-1] for c in cols.values())
    k = np.linspace(0, 2 * np.pi * fmax * b / U, nk)
    Q = np.full((nk, n, n), np.nan, complex)
    for j, (f, Qj) in cols.items():
        kj = 2 * np.pi * f * b / U
        for i in range(n):
            Q[:, i, j] = np.interp(k, kj, Qj[:, i].real) + 1j * np.interp(k, kj, Qj[:, i].imag)
    return {"mach": flow["Minf"], "U": U, "q": 0.5 * flow["rhoInf"] * U ** 2, "b": b, "k": k, "Q": Q,
            "modes": np.array(sorted(cols)), "missing": np.array(missing, int), "runs": np.array([str(r) for r in runs])}


def linearity(run_a, run_b, zero=None):
    """Relative difference of the Q column of two runs of the same mode at different amplitude."""
    ja, fa, Qa, _ = column(run_a, zero)
    jb, fb, Qb, _ = column(run_b, zero)
    if ja != jb:
        raise ValueError("the two runs move different modes")
    Qb = np.array([np.interp(fa, fb, Qb[:, i].real) + 1j * np.interp(fa, fb, Qb[:, i].imag) for i in range(Qb.shape[1])]).T
    scale = np.abs(Qa).max()
    return ja, fa, np.abs(Qb - Qa).max(axis=1) / scale


def load(path):
    d = np.load(path, allow_pickle=False)
    return {k: d[k] for k in d.files}


def table(paths):
    """Sorted list of GAF files by Mach."""
    g = sorted((load(p) for p in paths), key=lambda d: float(d["mach"]))
    return g
