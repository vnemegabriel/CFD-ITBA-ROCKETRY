"""ISA atmosphere and flight-envelope loading."""
import numpy as np


def isa(h):
    h = np.asarray(h, float)
    T = np.where(h < 11000, 288.15 - 0.0065 * h, 216.65)
    p = np.where(h < 11000, 101325 * (T / 288.15) ** 5.25588, 22632.1 * np.exp(-(h - 11000) / 6341.62))
    return p / (287.053 * T), np.sqrt(1.4 * 287.053 * T)


def load(csv, mach_min=1.05):
    """CSV columns: time [s], altitude [ft], vertical velocity [m/s]. Keeps points with M >= mach_min."""
    t, h_ft, v = np.loadtxt(csv, delimiter=",", comments="#").T
    rho, a = isa(h_ft * 0.3048)
    U = np.abs(v)
    k = U / a >= mach_min
    return {"t": t[k], "h": h_ft[k] * 0.3048, "rho": rho[k], "U": U[k], "mach": (U / a)[k], "q": 0.5 * rho[k] * U[k] ** 2}
