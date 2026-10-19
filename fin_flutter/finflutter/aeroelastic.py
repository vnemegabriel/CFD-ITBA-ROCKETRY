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


def roger_fit(k, Q, gamma):
    """Real matrices of Q(ik) ~ A0 + ik A1 + (ik)^2 A2 + sum_l ik/(ik + gamma_l) L_l, least squares over k."""
    k = np.asarray(k, float)
    s = 1j * k
    basis = np.column_stack([np.ones_like(s), s, s ** 2] + [s / (s + g) for g in gamma])
    X = np.vstack([basis.real, basis.imag])
    n = Q.shape[1]
    Y = np.vstack([Q.reshape(len(k), -1).real, Q.reshape(len(k), -1).imag])
    c = np.linalg.lstsq(X, Y, rcond=None)[0].reshape(-1, n, n)
    fit = np.einsum("km,mij->kij", basis, c)
    return {"A0": c[0], "A1": c[1], "A2": c[2], "lags": c[3:], "gamma": np.asarray(gamma, float),
            "residual": np.abs(fit - Q).max() / np.abs(Q).max()}


class RogerTable:
    """Roger fits at each CFD Mach number, interpolated linearly in Mach. Same gamma and b at every Mach."""

    def __init__(self, gafs, gamma=None):
        self.b = float(gafs[0]["b"])
        kmax = min(float(g["k"][-1]) for g in gafs)
        self.gamma = np.asarray(gamma if gamma is not None else kmax * np.array([0.2, 0.5, 1.0]))
        self.machs = np.array([float(g["mach"]) for g in gafs])
        self.fits = [roger_fit(g["k"], g["Q"], self.gamma) for g in gafs]

    def __call__(self, mach):
        if not self.machs[0] <= mach <= self.machs[-1]:
            return None
        i = min(np.searchsorted(self.machs, mach, side="right") - 1, len(self.machs) - 2)
        w = (mach - self.machs[i]) / (self.machs[i + 1] - self.machs[i]) if len(self.machs) > 1 else 0.0
        f0, f1 = self.fits[i], self.fits[min(i + 1, len(self.fits) - 1)]
        out = {key: (1 - w) * f0[key] + w * f1[key] for key in ("A0", "A1", "A2", "lags")}
        out.update(gamma=self.gamma, b=self.b)
        return out


def state_matrix(omega, g, rfa, q, U):
    """First-order form of the Roger system; states eta, eta', and one lag state vector per gamma."""
    n, b, L = len(omega), rfa["b"], len(rfa["gamma"])
    Mi = np.linalg.inv(np.eye(n) + q * (b / U) ** 2 * rfa["A2"])
    C = np.diag(g * omega) + q * (b / U) * rfa["A1"]
    K = np.diag(omega ** 2) + q * rfa["A0"]
    A = np.zeros((n * (2 + L), n * (2 + L)))
    A[:n, n:2 * n] = np.eye(n)
    A[n:2 * n, :n] = -Mi @ K
    A[n:2 * n, n:2 * n] = -Mi @ C
    for l in range(L):
        r = slice((2 + l) * n, (3 + l) * n)
        A[n:2 * n, r] = -q * Mi @ rfa["lags"][l]
        A[r, n:2 * n] = np.eye(n)
        A[r, r] = -(U / b) * rfa["gamma"][l] * np.eye(n)
    return A


def critical_q_fn(eig, scale, q_max, n_grid=200, tol=1e-6):
    """critical_q for any eig(q) -> eigenvalues. scale: a reference frequency [rad/s]."""
    unstable = lambda q: eig(q).real.max() > tol * scale
    qs = np.linspace(0, q_max, n_grid)
    hit = next((i for i, q in enumerate(qs) if unstable(q)), None)
    if hit is None:
        return np.inf, None, np.nan
    if hit == 0:
        return 0.0, "unstable at q = 0", np.nan
    lo, hi = qs[hit - 1], qs[hit]
    for _ in range(40):
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if unstable(mid) else (mid, hi)
    p = eig(hi)
    p = p[np.argmax(p.real)]
    return hi, ("divergence" if abs(p.imag) < 1e-3 * scale else "flutter"), abs(p.imag) / (2 * np.pi)


def envelope_rfa(omega, g, table, flight, q_factor=20.0):
    """envelope() with a RogerTable; flight points outside its Mach range are left out."""
    out = []
    for M, U, q in zip(flight["mach"], flight["U"], flight["q"]):
        rfa = table(M)
        if rfa is None:
            continue
        qc, kind, f = critical_q_fn(lambda qq: np.linalg.eigvals(state_matrix(omega, g, rfa, qq, U)), omega[0], q_factor * q)
        out.append({"mach": M, "U": U, "q": q, "q_crit": qc, "margin": qc / q, "kind": kind, "f_Hz": f})
    return out


def least_stable(omega, g, rfa, q, U):
    """Oscillatory root with the largest real part at (q, U), and its modal vector normalised to max |v| = 1, real there."""
    w, V = np.linalg.eig(state_matrix(omega, g, rfa, q, U))
    osc = np.flatnonzero(w.imag > 1e-3 * omega[0])
    i = osc[np.argmax(w.real[osc])]
    v = V[:len(omega), i]
    v = v / v[np.argmax(np.abs(v))]
    return w[i], v
