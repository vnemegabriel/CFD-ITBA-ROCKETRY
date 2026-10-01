"""Classical laminate theory for the T700/epoxy DB300 (+-45 NCF) / GA90R (woven 0/90) fin layups."""
import itertools
import numpy as np

EF1, EF2, GF12, NUF, RHO_F = 230e9, 15e9, 27e9, 0.20, 1800.0
EM, NUM, RHO_M = 3.5e9, 0.35, 1200.0
FAW = {"B": 0.300 / 2, "W": 0.302}  # kg/m2 per ply (DB300 counted per +-45 half)
KC = 0.92  # woven crimp knockdown

LAYUPS = {
    "ar1": "B+45 B-45 W0 W0 B+45 B-45 W0 B+45 B-45 W0 W0",
    "more_db": "B+45 B-45 B+45 B-45 W0 W0 B+45 B-45 W0 B+45 B-45 W0 B+45 B-45",
    "more_ga": "W0 W0 B+45 B-45 W0 W0 B+45 B-45 W0 W0 B+45 B-45 W0 W0",
    "equal": "B+45 B-45 W0 W0 B+45 B-45 W0 W0 B+45 B-45 W0 W0",
}


def _halpin_tsai(ef, em, vf, xi):
    eta = (ef / em - 1) / (ef / em + xi)
    return em * (1 + xi * eta * vf) / (1 - eta * vf)


def ply_props(vf):
    """(E1, E2, G12, nu12, t) per fabric at fibre volume fraction vf."""
    e1 = EF1 * vf + EM * (1 - vf)
    nu = NUF * vf + NUM * (1 - vf)
    e2 = _halpin_tsai(EF2, EM, vf, 2.0)
    g12 = _halpin_tsai(GF12, EM / (2 * (1 + NUM)), vf, 1.0)
    ew = 0.5 * (e1 + e2) * KC
    return {"B": (e1, e2, g12, nu, FAW["B"] / (vf * RHO_F)),
            "W": (ew, ew, g12 * KC, 2 * nu * e2 / (e1 + e2), FAW["W"] / (vf * RHO_F))}


def _qbar(e1, e2, g12, nu12, deg):
    d = 1 - nu12 ** 2 * e2 / e1
    q11, q22, q12, q66 = e1 / d, e2 / d, nu12 * e2 / d, g12
    m, n = np.cos(np.radians(deg)), np.sin(np.radians(deg))
    T = np.array([[m * m, n * n, 2 * m * n], [n * n, m * m, -2 * m * n], [-m * n, m * n, m * m - n * n]])
    Q = np.array([[q11, q12, 0], [q12, q22, 0], [0, 0, 2 * q66]])
    Qb = np.linalg.inv(T) @ Q @ T
    Qb[:, 2] /= 2
    return Qb


def sized_half_stack(layup, thickness_mm, vf):
    """Repeat the half-stack pattern (keeping +-45 pairs together) to the ply count closest to the mold thickness."""
    plies = [(s[0], float(s[1:])) for s in LAYUPS[layup].split()]
    groups, i = [], 0
    while i < len(plies):
        n = 2 if plies[i][0] == "B" and i + 1 < len(plies) and plies[i + 1][0] == "B" else 1
        groups.append(plies[i:i + n]); i += n
    props, half, h, target = ply_props(vf), [], 0.0, thickness_mm * 5e-4
    for g in itertools.cycle(groups):
        dt = sum(props[f][4] for f, _ in g)
        if h + dt > target:
            if abs(target - h - dt) < abs(target - h) or not half:
                half += g
            break
        half += g; h += dt
    return half


def laminate(layup="ar1", thickness_mm=8.0, vf=0.5, beta_deg=0.0):
    """Closed-mould laminate: true Vf from areal weight, plies rotated by beta. Returns D (3x3, x = chord), t, rho, G13."""
    half = sized_half_stack(layup, thickness_mm, vf)
    t = thickness_mm * 1e-3
    vf = 2 * sum(FAW[f] for f, _ in half) / (RHO_F * t)
    props = ply_props(vf)
    stack = [(f, a + beta_deg) for f, a in half + half[::-1]]
    z, D, G = -t / 2, np.zeros((3, 3)), 0.0
    for f, a in stack:
        e1, e2, g12, nu, tp = props[f]
        D += _qbar(e1, e2, g12, nu, a) * ((z + tp) ** 3 - z ** 3) / 3
        G += g12 * tp; z += tp
    return {"D": D, "t": t, "rho": vf * RHO_F + (1 - vf) * RHO_M, "G13": G / t, "vf": vf}


def isotropic(E, nu, t, rho):
    D = E * t ** 3 / (12 * (1 - nu ** 2)) * np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]])
    return {"D": D, "t": t, "rho": rho, "G13": E / (2 * (1 + nu))}
