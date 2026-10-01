"""MITC4 Mindlin plate FE (DOFs per node: w, bx, by with u = z*bx, v = z*by) and modal analysis."""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

_G = np.array([-1, 1]) / np.sqrt(3)
_XI = np.array([-1, 1, 1, -1]), np.array([-1, -1, 1, 1])


def mesh(cr, ct, span, sweep_deg, nx, ny):
    """Structured Q4 mesh of a trapezoidal fin: x chordwise (flow), y spanwise, root at y = 0."""
    s, e = np.meshgrid(np.linspace(0, 1, nx + 1), np.linspace(0, 1, ny + 1))
    y = e * span
    x = y * np.tan(np.radians(sweep_deg)) + s * (cr + (ct - cr) * e)
    n = np.arange((nx + 1) * (ny + 1)).reshape(ny + 1, nx + 1)
    elems = np.stack([n[:-1, :-1], n[:-1, 1:], n[1:, 1:], n[1:, :-1]], -1).reshape(-1, 4)
    return np.column_stack([x.ravel(), y.ravel()]), elems


def _shape(xi, eta):
    N = (1 + _XI[0] * xi) * (1 + _XI[1] * eta) / 4
    dN = np.array([_XI[0] * (1 + _XI[1] * eta), _XI[1] * (1 + _XI[0] * xi)]) / 4
    return N, dN


def _shear_nat(xy, xi, eta):
    """Rows: covariant transverse shear strains (g_xi, g_eta) at a point, as a 2x12 operator."""
    N, dN = _shape(xi, eta)
    J = dN @ xy
    B = np.zeros((2, 12))
    B[:, 0::3] = dN
    B[:, 1::3] = J[:, :1] * N
    B[:, 2::3] = J[:, 1:] * N
    return B


def _element(xy, D, Hs, rho, t):
    K, M, Sx, Mw = np.zeros((12, 12)), np.zeros((12, 12)), np.zeros((4, 4)), np.zeros((4, 4))
    tie = {k: _shear_nat(xy, *p) for k, p in {"A": (0, 1), "C": (0, -1), "B": (-1, 0), "D": (1, 0)}.items()}
    for xi in _G:
        for eta in _G:
            N, dN = _shape(xi, eta)
            J = dN @ xy
            dA = np.linalg.det(J)
            dx = np.linalg.solve(J, dN)
            Bb = np.zeros((3, 12))
            Bb[0, 1::3], Bb[1, 2::3] = dx[0], dx[1]
            Bb[2, 1::3], Bb[2, 2::3] = dx[1], dx[0]
            gnat = np.vstack([(1 + eta) / 2 * tie["A"][0] + (1 - eta) / 2 * tie["C"][0],
                              (1 + xi) / 2 * tie["D"][1] + (1 - xi) / 2 * tie["B"][1]])
            Bs = np.linalg.solve(J, gnat)
            K += (Bb.T @ D @ Bb + Bs.T @ Hs @ Bs) * dA
            Nw = np.zeros((3, 12))
            for i in range(3):
                Nw[i, i::3] = N
            M += Nw.T @ np.diag([rho * t, rho * t ** 3 / 12, rho * t ** 3 / 12]) @ Nw * dA
            Sx += np.outer(N, dx[0]) * dA
            Mw += np.outer(N, N) * dA
    return K, M, Sx, Mw


def assemble(nodes, elems, lam, kappa=5 / 6):
    """Global K, M (3 DOF/node) and w-DOF aero operators Sx = int N_i dN_j/dx, Mw = int N_i N_j."""
    Hs = kappa * lam["G13"] * lam["t"] * np.eye(2)
    out = {k: ([], [], []) for k in "KMSW"}
    for e in elems:
        dof = (3 * e[:, None] + np.arange(3)).ravel()
        for k, m, d in zip("KMSW", _element(nodes[e], lam["D"], Hs, lam["rho"], lam["t"]), (dof, dof, e, e)):
            r, c, v = out[k]
            r.append(np.repeat(d, len(d))); c.append(np.tile(d, len(d))); v.append(m.ravel())
    size = {"K": 3 * len(nodes), "M": 3 * len(nodes), "S": len(nodes), "W": len(nodes)}
    return {k: sp.csr_matrix((np.concatenate(v), (np.concatenate(r), np.concatenate(c))), shape=(size[k],) * 2)
            for k, (r, c, v) in out.items()}


def modes(mats, fixed_dofs, n):
    """Lowest n modes. Returns omega [rad/s] and mass-normalised w-components phi_w [nNodes x n]."""
    free = np.setdiff1d(np.arange(mats["K"].shape[0]), fixed_dofs)
    K, M = mats["K"][free][:, free], mats["M"][free][:, free]
    w2, v = eigsh(K, n, M, sigma=-1.0)
    i = np.argsort(w2)
    v = v[:, i] / np.sqrt(np.einsum("ij,ij->j", v[:, i], M @ v[:, i]))
    phi = np.zeros((mats["K"].shape[0], n))
    phi[free] = v
    return np.sqrt(np.abs(w2[i])), phi[0::3]


def clamped_dofs(nodes, x_ranges, tol=1e-9):
    """All DOFs of root nodes (y = 0) whose x lies in any [x0, x1] interval."""
    x, y = nodes[:, 0], nodes[:, 1]
    sel = np.flatnonzero((y < tol) & np.any([(x >= a - tol) & (x <= b + tol) for a, b in x_ranges], axis=0))
    return (3 * sel[:, None] + np.arange(3)).ravel()
