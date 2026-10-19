"""Figures for docs/PROCEDURE.md. Each function writes one PNG and returns its path."""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3, "figure.dpi": 110})


def _save(fig, out):
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", dpi=130)
    plt.close(fig)
    return out


def geometry(fin, out):
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.6), gridspec_kw={"width_ratios": [1.1, 1]})
    a.plot([0, fin.cr, fin.sweep + fin.ct, fin.sweep, 0], [0, 0, fin.span, fin.span, 0], "k-")
    for e in (0.0, 0.5, 1.0):
        x0, c = fin.x_le(e), fin.chord(e)
        a.plot([x0, x0 + c], [e * fin.span] * 2, "--", lw=0.8, label=f"eta = {e:.1f}")
    a.set_aspect("equal"); a.set_xlabel("x desde el LE de raíz [m]"); a.set_ylabel("y desde la raíz [m]")
    a.set_title(f"Planta: cr {fin.cr:.3f}, ct {fin.ct:.3f}, span {fin.span:.3f} m, flecha LE {fin.sweep_deg:.1f}°")
    a.legend(loc="lower right", fontsize=7)
    eta, xi = np.meshgrid(np.linspace(0, 1, 81), np.linspace(0, 1, 401))
    x, y = fin.x_le(eta) + xi * fin.chord(eta), eta * fin.span
    t = np.vectorize(fin.thickness)(x, y) * 1e3
    cs = b.contourf(x, y, t, levels=np.arange(0, fin.t * 1e3 + 0.5, 1.0), cmap="viridis")
    fig.colorbar(cs, ax=b, label="espesor local [mm]")
    b.set_aspect("equal"); b.set_xlabel("x [m]"); b.set_ylabel("y [m]")
    b.set_title(f"Espesor: {fin.t*1e3:.1f} mm a c/2, {2*fin.rho*1e3:.0f} mm en LE/TE")
    return _save(fig, out)


def modes(nodes, elems, phi, omega, out, n=4, kinds=None):
    tri = mtri.Triangulation(nodes[:, 0], nodes[:, 1], np.vstack([elems[:, [0, 1, 2]], elems[:, [0, 2, 3]]]))
    fig, axs = plt.subplots(1, n, figsize=(3.0 * n, 3.0), sharey=True)
    for k, ax in enumerate(axs):
        w = phi[:, k] / np.abs(phi[:, k]).max()
        cs = ax.tricontourf(tri, w, levels=np.linspace(-1, 1, 21), cmap="RdBu_r")
        ax.tricontour(tri, w, levels=[0], colors="k", linewidths=0.8)
        ax.set_aspect("equal"); ax.set_title(f"modo {k+1}: {omega[k]/2/np.pi:.0f} Hz" + (f"\n{kinds[k]}" if kinds else ""))
        ax.set_xlabel("x [m]")
    axs[0].set_ylabel("y [m]")
    fig.colorbar(cs, ax=axs, shrink=0.8, label="w / max|w|  (línea negra: nodal)")
    return _save(fig, out)


def convergence(rows, out):
    """rows: [(label, nElems, [f1..fn])]"""
    fig, ax = plt.subplots(figsize=(5, 3.2))
    f = np.array([r[2] for r in rows]); ne = [r[1] for r in rows]
    for k in range(f.shape[1]):
        ax.plot(ne, f[:, k] / f[-1, k], "o-", label=f"modo {k+1}")
    ax.set_xscale("log"); ax.set_xlabel("elementos"); ax.set_ylabel("f / f(malla más fina)")
    ax.set_title("Convergencia de frecuencias con la malla FE"); ax.legend(fontsize=7)
    return _save(fig, out)


def gaf(d, out, n=4):
    """Re and Im of Q_ij(k), i, j <= n."""
    k, Q = d["k"], d["Q"]
    n = min(n, Q.shape[1])
    fig, axs = plt.subplots(n, n, figsize=(2.6 * n, 2.2 * n), sharex=True)
    for i in range(n):
        for j in range(n):
            ax = axs[i, j]
            ax.plot(k, Q[:, i, j].real, label="Re")
            ax.plot(k, Q[:, i, j].imag, label="Im")
            ax.set_title(f"Q{i+1}{j+1}", fontsize=8)
            if i == n - 1:
                ax.set_xlabel("k = ωb/U")
    axs[0, 0].legend(fontsize=7)
    fig.suptitle(f"GAF, M = {float(d['mach']):.3f}  (b = {float(d['b']):.4f} m)")
    return _save(fig, out)


def linearity(f, err, j, out):
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.plot(f, err * 100)
    ax.set_xlabel("f [Hz]"); ax.set_ylabel("|Q(2A) − Q(A)| / max|Q| [%]")
    ax.set_title(f"Linealidad, columna del modo {j+1}")
    return _save(fig, out)


def envelope(env, bound, out):
    fig, ax = plt.subplots(figsize=(6, 3.4))
    m = np.array([r["mach"] for r in env]); mg = np.array([r["margin"] for r in env])
    mb = np.array([r["mach"] for r in bound]); gb = np.array([r["margin"] for r in bound])
    ax.semilogy(mb, np.where(np.isfinite(gb), gb, np.nan), ".", ms=2, color="0.6", label="piston (Ackeret)")
    ax.semilogy(m, np.where(np.isfinite(mg), mg, np.nan), ".", ms=3, label="CFD + Roger")
    ax.axhline(1.0, color="r", lw=0.8); ax.axhline(1.15 ** 2, color="r", lw=0.8, ls="--", label="15 % en velocidad")
    ax.set_xlabel("M"); ax.set_ylabel("q_crit / q (a Mach constante)"); ax.legend(fontsize=7)
    return _save(fig, out)


def tap(rows, f_model, out):
    fig, (a, b) = plt.subplots(1, 2, figsize=(9, 3.2))
    f = np.array([r["f_Hz"] for r in rows]); g = np.array([r["g"] for r in rows]); amp = np.array([r["amplitude"] for r in rows])
    a.plot(f, g, "o", ms=3)
    for fm in f_model:
        a.axvline(fm, color="0.6", lw=0.8)
    a.set_xlabel("f [Hz] (gris: modelo FE)"); a.set_ylabel("g = 2ζ")
    sc = b.scatter(amp, g, c=f, s=10, cmap="viridis")
    fig.colorbar(sc, ax=b, label="f [Hz]")
    b.set_xlabel("amplitud del golpe [u.a.]"); b.set_ylabel("g"); b.set_title("¿g cambia con la amplitud?")
    return _save(fig, out)


def lco(rows, cross, out):
    fig, ax = plt.subplots(figsize=(5, 3.2))
    A = np.array([r["A"] for r in rows]) * 1e3
    ax.plot(A, [r["W_a"] for r in rows], "o-", label="W_a (aire → aleta)")
    ax.plot(A, [r["W_d"] for r in rows], "s-", label="W_d (disipación estructural)")
    for x in cross:
        ax.axvline(x["A"] * 1e3, color="r", ls="--", lw=0.8)
    ax.set_xlabel("A [mm]"); ax.set_ylabel("energía por ciclo [J]"); ax.legend(fontsize=7)
    return _save(fig, out)
