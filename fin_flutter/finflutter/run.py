import json
import subprocess
from pathlib import Path

import numpy as np
import yaml

from . import aero, aeroelastic, flight, geometry, laminate, structure


def load_case(case):
    return yaml.safe_load(Path(case).read_text()) if not isinstance(case, dict) else case


def build_structure(c):
    """Fin, FE mesh, matrices and modes of a case. Returns a dict."""
    gm, lm = c["geometry"], dict(c["laminate"])
    thickness = None
    if gm.get("source") == "flightModel":
        fin = geometry.fin_from_flightmodel(dims=gm.get("dims"))
        planform = dict(cr=fin.cr, ct=fin.ct, span=fin.span, sweep_deg=fin.sweep_deg)
        lm["thickness_mm"] = fin.t * 1e3
        if lm.pop("taper", "scaled") == "scaled":
            thickness = fin.thickness
    else:
        fin, planform = None, gm
    lam = laminate.laminate(**lm)
    nodes, elems = structure.mesh(**planform, **c["mesh"])
    mats = structure.assemble(nodes, elems, lam, thickness=thickness)
    clamp = [[0.0, planform["cr"]]] if c.get("clamp_x", "root") == "root" else c["clamp_x"]
    omega, phi = structure.modes(mats, structure.clamped_dofs(nodes, clamp), c["modes"])
    return {"fin": fin, "planform": planform, "lam": lam, "nodes": nodes, "elems": elems, "mats": mats,
            "omega": omega, "phi": phi, "thickness": thickness}


def run(case, root="."):
    c = load_case(case)
    s = build_structure(c)
    omega, phi, mats = s["omega"], s["phi"], s["mats"]
    pts = flight.load(Path(root) / c["flight"]["csv"], c["flight"].get("mach_min", 1.05))
    env = aeroelastic.envelope(omega, c["g"] * np.ones_like(omega), aero.PistonTheory(phi, mats, **c["aero"]), pts,
                               c.get("q_factor", 20.0))
    crit = min(env, key=lambda r: r["margin"])
    return {"case": c, "f_Hz": (omega / 2 / np.pi).tolist(), "D": s["lam"]["D"].tolist(), "critical": crit,
            "envelope": env}


def git_sha():
    return subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()


def main(case):
    case = Path(case)
    r = run(case)
    c = r["critical"]
    print("f [Hz]:", " ".join(f"{f:.0f}" for f in r["f_Hz"]))
    print(f"min margin q_crit/q = {c['margin']:.2f} ({c['kind'] or 'stable'}, f={c['f_Hz']:.0f} Hz) at M={c['mach']:.3f}, "
          f"q={c['q']:.0f} Pa; constant-Mach velocity margin {np.sqrt(c['margin']):.2f}")
    out = Path("results") / f"{case.stem}.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps({"commit": git_sha(), **r}, indent=1, default=float))
    print("->", out)
