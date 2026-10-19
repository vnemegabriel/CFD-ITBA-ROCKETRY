#!/usr/bin/env python3
"""Cd of the IREC flight from the TeleMega, against the CFD runs and OpenRocket.

    cdFlight.py [--csv AconcaguaIREC_TeleMega.csv] [--ork Aconcagua.ork] [--mass <kg>]
                [--runs ~/runs] [--out cdFlight.png]

The CSV uses a dot as decimal sign and a comma between fields, with the column
names of the TeleMega export (AltOS); load() stops on a cell outside the
physical limits of LIMITS.  The mass after burnout comes from the OpenRocket
simulation in the .ork file, or from --mass.

Everything else comes from the telemetry:

- Atmosphere of Midland: pad pressure p0 is the mean before launch.  T0 and the
  lapse rate come from fitting p = p0 (1 - L h/T0)^(g/RL) to the GPS height
  above the pad (sqrt(pad_range^2 - pad_dist^2)) against the pressure, on the
  descent, where the baro sees no shock.
- Speed along the path: dV/dt = accel_x - g cos(tilt), integrated from an
  anchor in the subsonic coast where V = (dh/dt from the baro) / cos(tilt).
  The baro and AltOS speed are not used elsewhere: past M 0.8 the static port
  sits behind shocks, and AltOS converts pressure with the ISA.
- Height by integrating V cos(tilt) from the anchor; ambient p and T from that
  height and the fitted profile.
- Drag D = -m accel_x in the coast (accel_x is the axial specific force).
"""
import argparse
import csv
import glob
import math
import os
import re
import statistics
import xml.etree.ElementTree as ET
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

R, GAMMA, G = 287.05, 1.4, 9.7915
AREF = math.pi * 0.0755 ** 2
COLS = ("time", "state_name", "pressure", "speed", "accel_x", "tilt",
        "pad_range", "pad_dist", "nsat")


# Physical limits of a sounding-rocket flight; a value outside them is a unit,
# format or sensor error, not data.
LIMITS = {"time": (-60, 4000), "pressure": (1e4, 1.2e5), "speed": (-400, 1000),
          "accel_x": (-400, 400), "tilt": (0, 180), "pad_range": (0, 5e4),
          "pad_dist": (0, 5e4), "nsat": (0, 40)}


def num(s):
    s = s.strip()
    if "," in s:
        raise ValueError("decimal comma in %r: the CSV must use a dot as decimal sign" % s)
    return float(s) if s else float("nan")


def load(path):
    with open(path, encoding="utf-8") as f:
        r = csv.reader(f)
        head = next(r)
        head[0] = head[0].lstrip("#")
        idx = {k: head.index(k) for k in COLS}
        rows = []
        for line, row in enumerate(r, 2):
            if len(row) <= max(idx.values()):
                continue
            d = {k: (row[i] if k == "state_name" else num(row[i])) for k, i in idx.items()}
            for k, (lo, hi) in LIMITS.items():
                if not lo <= d[k] <= hi:
                    raise ValueError("%s line %d: %s = %g is outside [%g, %g]" % (path, line, k, d[k], lo, hi))
            rows.append(d)
    return rows


def atmosphere(rows):
    pre = [d["pressure"] for d in rows if d["time"] < 0]
    p0 = sum(pre) / len(pre)
    pts, prev = [], None
    for d in rows:
        if d["time"] < 30 or d["nsat"] < 6:
            continue
        k = (d["pad_range"], d["pad_dist"])
        if k != prev and d["pad_range"] > d["pad_dist"]:
            pts.append((math.sqrt(d["pad_range"] ** 2 - d["pad_dist"] ** 2), d["pressure"]))
        prev = k
    best = None
    for T0 in [x * 0.25 for x in range(1160, 1280)]:
        for L in [x * 0.0001 for x in range(40, 95)]:
            e = sum((math.log(p) - math.log(p0 * (1 - L * h / T0) ** (G / (R * L)))) ** 2
                    for h, p in pts)
            if best is None or e < best[0]:
                best = (e, T0, L)
    e, T0, L = best
    return p0, T0, L, math.sqrt(e / len(pts)), len(pts)


def trajectory(rows, atm, tanchor, tend):
    p0, T0, L = atm[:3]

    def h_of_p(p):
        return T0 / L * (1 - (p / p0) ** (R * L / G))

    fl = [d for d in rows if 3.0 <= d["time"] <= tend]
    win = [d for d in fl if abs(d["time"] - tanchor) <= 0.5]
    a = [d for d in win if d["time"] < tanchor - 0.3]
    b = [d for d in win if d["time"] > tanchor + 0.3]
    ta, tb = (sum(d["time"] for d in s) / len(s) for s in (a, b))
    ha, hb = (h_of_p(sum(d["pressure"] for d in s) / len(s)) for s in (a, b))
    i0 = min(range(len(fl)), key=lambda i: abs(fl[i]["time"] - tanchor))
    w0 = (hb - ha) / (tb - ta)
    V = [0.0] * len(fl)
    H = [0.0] * len(fl)
    V[i0] = w0 / math.cos(math.radians(fl[i0]["tilt"]))
    H[i0] = h_of_p(fl[i0]["pressure"])

    def dVdt(d):
        return d["accel_x"] - G * math.cos(math.radians(d["tilt"]))

    for i in range(i0 + 1, len(fl)):
        dt = fl[i]["time"] - fl[i - 1]["time"]
        V[i] = V[i - 1] + 0.5 * (dVdt(fl[i]) + dVdt(fl[i - 1])) * dt
    for i in range(i0 - 1, -1, -1):
        dt = fl[i + 1]["time"] - fl[i]["time"]
        V[i] = V[i + 1] - 0.5 * (dVdt(fl[i]) + dVdt(fl[i + 1])) * dt
    w = [v * math.cos(math.radians(d["tilt"])) for v, d in zip(V, fl)]
    for i in range(i0 + 1, len(fl)):
        H[i] = H[i - 1] + 0.5 * (w[i] + w[i - 1]) * (fl[i]["time"] - fl[i - 1]["time"])
    for i in range(i0 - 1, -1, -1):
        H[i] = H[i + 1] - 0.5 * (w[i] + w[i + 1]) * (fl[i + 1]["time"] - fl[i]["time"])
    out = []
    for d, v, h in zip(fl, V, H):
        T = T0 - L * h
        p = p0 * (T / T0) ** (G / (R * L))
        out.append(dict(t=d["time"], state=d["state_name"], V=v, h=h, T=T, p=p,
                        M=v / math.sqrt(GAMMA * R * T), ax=d["accel_x"], tilt=d["tilt"]))
    return out


def drag(traj, mass, tstart, vmin):
    out = []
    for s in traj:
        if s["t"] < tstart or s["V"] < vmin or s["ax"] >= 0 or s["tilt"] >= 80:
            continue
        rho = s["p"] / (R * s["T"])
        out.append((s["t"], s["M"], -mass * s["ax"] / (0.5 * rho * s["V"] ** 2 * AREF)))
    return out


def binned(pts, width):
    bins = {}
    for _, m, cd in pts:
        bins.setdefault(round(m / width), []).append(cd)
    res = []
    for k in sorted(bins):
        v = sorted(bins[k])
        if len(v) >= 5:
            q = statistics.quantiles(v, n=4)
            res.append((k * width, statistics.median(v), q[0], q[2], len(v)))
    return res


def cfd(runs):
    out = []
    for d in sorted(glob.glob(os.path.join(runs, "*-a[0-9][0-9]-*"))):
        info = os.path.join(d, "run.info")
        coef = glob.glob(os.path.join(d, "postProcessing/forceCoeffs1/*/coefficient.dat"))
        if not os.path.exists(info) or not coef:
            continue
        txt = open(info).read()
        opts = re.search(r'options\s+"([^"]*)"', txt).group(1)
        minf = float(re.search(r"--Minf (\S+)", opts).group(1))
        mesh = re.search(r"--mesh (\S+)", opts).group(1).rstrip("/").split("/")[-1]
        solver = re.search(r"--solver (\S+)", opts)
        solver = solver.group(1) if solver else "central"
        alpha = re.search(r"--alpha (\S+)", opts)
        alpha = float(alpha.group(1)) if alpha else 0.0
        cd = [float(l.split()[1]) for l in open(max(coef, key=lambda p: float(p.split(os.sep)[-2]))) if l.strip() and not l.startswith("#")]
        if len(cd) < 20 or "status      ok" not in txt:
            continue
        out.append((os.path.basename(d), minf, alpha, mesh, solver, statistics.mean(cd[-20:])))
    return out


def openrocket(path):
    with zipfile.ZipFile(path) as z:
        name = [n for n in z.namelist() if n.endswith(".ork")][0]
        root = ET.fromstring(z.read(name))
    sim = next(root.iter("simulation"))
    burnout = min(float(e.get("time")) for e in sim.iter("event") if e.get("type") == "burnout")
    apogee = min(float(e.get("time")) for e in sim.iter("event") if e.get("type") == "apogee")
    br = next(sim.iter("databranch"))
    types = br.get("types").split(",")
    it, im, ic = types.index("Time"), types.index("Mach number"), types.index("Drag coefficient")
    imass = types.index("Mass")
    pts, mass = [], None
    for p in br.iter("datapoint"):
        v = p.text.split(",")
        if burnout < float(v[it]) < apogee:
            pts.append((float(v[im]), float(v[ic])))
            if mass is None:
                mass = float(v[imass])
    return sorted(pts), mass


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=os.path.join(here, "AconcaguaIREC_TeleMega.csv"))
    ap.add_argument("--mass", type=float, help="mass after burnout [kg]; default: from --ork")
    ap.add_argument("--runs", default=os.path.expanduser("~/runs"))
    ap.add_argument("--ork", default=os.path.join(here, "Aconcagua.ork"))
    ap.add_argument("--vmin", type=float, default=60)
    ap.add_argument("--tstart", type=float, default=4.0)
    ap.add_argument("--tanchor", type=float, default=18.0)
    ap.add_argument("--alphas", type=lambda s: [float(x) for x in s.split(",")], default=[0.0, 8.0])
    ap.add_argument("--out", default=os.path.join(here, "cdFlight.png"))
    a = ap.parse_args()
    ork, ork_mass = openrocket(a.ork) if os.path.exists(a.ork) else ([], None)
    if a.mass is None:
        if ork_mass is None:
            ap.error("give --mass, or an --ork file with a simulation")
        a.mass = ork_mass
    print("mass after burnout: %.2f kg" % a.mass)

    rows = load(a.csv)
    atm = atmosphere(rows)
    p0, T0, L, rms, n = atm
    print("Midland: p0 %.0f Pa (pad, before launch), T0 %.2f K, lapse %.2f K/km "
          "(%d GPS points, rms ln p %.4f)" % (p0, T0, L * 1000, n, rms))
    gps_apo = max((math.sqrt(d["pad_range"] ** 2 - d["pad_dist"] ** 2), d["time"])
                  for d in rows if 30 < d["time"] < 45 and d["nsat"] >= 6
                  and d["pad_range"] > d["pad_dist"])
    traj = trajectory(rows, atm, a.tanchor, gps_apo[1])
    top = max(traj, key=lambda s: s["h"])
    print("check: integrated apogee %.0f m at t %.2f;  GPS %.0f m at t %.2f" % (
        top["h"], top["t"], gps_apo[0], gps_apo[1]))

    print("\nsweep conditions in the coast (Midland atmosphere):")
    print("%5s %7s %8s %9s %8s" % ("M", "t [s]", "h pad", "pInf", "Tinf"))
    for m in (1.8, 1.2, 0.8, 0.6, 0.3):
        s = min((s for s in traj if s["t"] >= a.tstart), key=lambda s: abs(s["M"] - m))
        print("%5.1f %7.2f %8.0f %9.0f %8.2f" % (m, s["t"], s["h"], s["p"], s["T"]))

    tel = drag(traj, a.mass, a.tstart, a.vmin)
    bins = binned(tel, 0.025)
    runs = [r for r in cfd(a.runs) if r[2] in a.alphas]

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.scatter([m for _, m, _ in tel], [c for _, _, c in tel], s=3, color="0.75",
               label="TeleMega, each sample")
    ax.plot([b[0] for b in bins], [b[1] for b in bins], color="k", lw=2,
            label="TeleMega, median in bins of ΔM = 0.025")
    ax.fill_between([b[0] for b in bins], [b[2] for b in bins], [b[3] for b in bins],
                    color="k", alpha=0.15, lw=0, label="TeleMega, quartiles")
    if ork:
        ax.plot([m for m, _ in ork], [c for _, c in ork], color="tab:green", lw=1.5,
                ls="--", label="OpenRocket")
    marker = {"h0": "o", "h1": "s", "h2": "D", "h3": "^"}
    color = {0: "tab:blue", 8: "tab:orange", 16: "tab:red"}
    shift = {"h0": -0.012, "h1": 0.0, "h2": 0.012, "h3": 0.024}
    for name, m, al, mesh, solver, cd in runs:
        hollow = m > 1 and solver == "pimple"
        col = color.get(int(al), "k")
        ax.plot(m + shift.get(mesh, 0), cd, marker.get(mesh, "x"), ms=8,
                mfc="none" if hollow else col, mec=col, mew=2, ls="none")
    handles, labels = ax.get_legend_handles_labels()
    for al in sorted({int(r[2]) for r in runs}):
        handles.append(Line2D([], [], color=color.get(al, "k"), marker="o", ls="none", ms=8))
        labels.append("CFD α = %d°" % al)
    for mesh in sorted({r[3] for r in runs}):
        handles.append(Line2D([], [], color="0.3", marker=marker.get(mesh, "x"), ls="none",
                              mfc="none", ms=8))
        labels.append("mesh %s" % mesh)
    if any(r[1] > 1 and r[4] == "pimple" for r in runs):
        handles.append(Line2D([], [], color="0.3", marker="o", ls="none", mfc="none", mew=2, ms=8))
        labels.append("hollow: supersonic pimple (start run only)")
    ax.set_xlabel("Mach")
    ax.set_ylabel("Cd  (Aref = π·0.0755² m²)")
    ax.set_title("Aconcagua, IREC Midland: flight Cd (coast) against CFD, α = %s°"
                 % ", ".join("%g" % x for x in a.alphas))
    ax.set_xlim(0.2, 1.95)
    ax.set_ylim(0, 0.9)
    ax.grid(alpha=0.3)
    ax.legend(handles, labels, fontsize=8, loc="upper right", ncol=2)
    ax.text(0.99, 0.02, "mass %.2f kg (OpenRocket); Midland atmosphere fitted to the GPS: "
            "T0 = %.1f K, %.1f K/km, p0 = %.0f Pa" % (a.mass, T0, L * 1000, p0),
            transform=ax.transAxes, ha="right", fontsize=7, color="0.4")
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)

    print("\n%-20s %5s %5s %4s %8s %8s %8s %7s" % (
        "case", "M", "alpha", "mesh", "CFD", "flight", "OR", "CFD/fl"))
    for name, m, al, mesh, solver, cd in runs:
        near = [b[1] for b in bins if abs(b[0] - m) < 0.0126]
        orv = min(ork, key=lambda p: abs(p[0] - m))[1] if ork else float("nan")
        fl = near[0] if near else float("nan")
        print("%-20s %5.2f %5.0f %4s %8.3f %8.3f %8.3f %7.2f" % (
            name, m, al, mesh, cd, fl, orv, cd / fl))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
