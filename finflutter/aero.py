"""Generalised aerodynamic force providers.

A provider is callable as provider(mach) -> (A0, B1), the modal aero force per unit dynamic pressure being
    f = -q * (A0 @ eta + B1 @ eta_dot / U).
"""
import numpy as np


class PistonTheory:
    """Linear supersonic quasi-steady theory on both fin faces, Delta p = sides * c0 * q * (w_x + c1/c0 * w_t / U).

    'ackeret': c0 = 2/beta, c1 = c0 (M^2-2)/(M^2-1)  (Van Dyke, Dowell);  'lighthill': c0 = c1 = 2/M.
    """

    def __init__(self, phi_w, mats, theory="ackeret", sides=2):
        self.S, self.W = phi_w.T @ (mats["S"] @ phi_w), phi_w.T @ (mats["W"] @ phi_w)
        self.theory, self.sides = theory, sides

    def __call__(self, mach):
        if self.theory == "lighthill":
            c0 = c1 = 2 / mach
        else:
            c0 = 2 / np.sqrt(mach ** 2 - 1)
            c1 = c0 * (mach ** 2 - 2) / (mach ** 2 - 1)
        return self.sides * c0 * self.S, self.sides * c1 * self.W


class TabulatedGAF:
    """GAFs tabulated from CFD/panel codes: Q[m, k] (n x n, per unit q) at reduced frequencies k = omega*b/U.

    Fits Q(ik) = A0 + ik*A1 per Mach in least squares, then interpolates linearly in Mach; B1 = b*A1.
    """

    def __init__(self, machs, ks, Q, b):
        V = np.column_stack([np.ones(len(ks)), 1j * np.asarray(ks)])
        Q = np.asarray(Q)
        coef = np.stack([np.linalg.lstsq(V, Qm.reshape(len(ks), -1), rcond=None)[0] for Qm in Q])
        n = Q.shape[-1]
        self.machs, self.A0, self.B1 = np.asarray(machs), coef[:, 0].real.reshape(-1, n, n), b * coef[:, 1].real.reshape(-1, n, n)

    def __call__(self, mach):
        f = lambda A: np.stack([np.interp(mach, self.machs, a) for a in A.reshape(len(A), -1).T]).reshape(A.shape[1:])
        return f(self.A0), f(self.B1)
