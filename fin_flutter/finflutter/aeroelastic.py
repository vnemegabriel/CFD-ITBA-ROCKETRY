"""Exact modal aeroelastic stability: eta'' + (g*Om + q*B1/U) eta' + (Om^2 + q*A0) eta = 0."""
import numpy as np


def roots(omega, g, A0, B1, q, U):
    n = len(omega)
    K = np.diag(omega ** 2) + q * A0
    C = np.diag(g * omega) + q * B1 / U
    return np.linalg.eigvals(np.block([[np.zeros((n, n)), np.eye(n)], [-K, -C]]))


def critical_q(omega, g, A0, B1, U, q_max, n_grid=200, tol=1e-6):
    """Lowest q in (0, q_max] with an unstable root. Returns (q, 'flutter'|'divergence', freq_Hz) or (inf, None, nan)."""
    unstable = lambda q: roots(omega, g, A0, B1, q, U).real.max() > tol * omega[0]
    qs = np.linspace(0, q_max, n_grid)
    hit = next((i for i, q in enumerate(qs) if unstable(q)), None)
    if hit is None:
        return np.inf, None, np.nan
    lo, hi = qs[hit - 1], qs[hit]
    for _ in range(40):
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if unstable(mid) else (mid, hi)
    p = roots(omega, g, A0, B1, hi, U)
    p = p[np.argmax(p.real)]
    return hi, ("divergence" if abs(p.imag) < 1e-3 * omega[0] else "flutter"), abs(p.imag) / (2 * np.pi)


def envelope(omega, g, aero, flight, q_factor=20.0):
    """Constant-Mach stability margin q_crit/q_flight at every flight point."""
    out = []
    for M, U, q in zip(flight["mach"], flight["U"], flight["q"]):
        A0, B1 = aero(M)
        qc, kind, f = critical_q(omega, g, A0, B1, U, q_factor * q)
        out.append({"mach": M, "U": U, "q": q, "q_crit": qc, "margin": qc / q, "kind": kind, "f_Hz": f})
    return out
