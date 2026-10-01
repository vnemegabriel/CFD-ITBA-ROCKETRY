import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml

from . import aero, aeroelastic, flight, laminate, structure


def run(case, root="."):
    c = yaml.safe_load(Path(case).read_text()) if not isinstance(case, dict) else case
    gm, lm = c["geometry"], c["laminate"]
    lam = laminate.laminate(**lm)
    nodes, elems = structure.mesh(**gm, **c["mesh"])
    mats = structure.assemble(nodes, elems, lam)
    omega, phi = structure.modes(mats, structure.clamped_dofs(nodes, c["clamp_x"]), c["modes"])
    pts = flight.load(Path(root) / c["flight"]["csv"], c["flight"].get("mach_min", 1.05))
    env = aeroelastic.envelope(omega, c["g"] * np.ones_like(omega), aero.PistonTheory(phi, mats, **c["aero"]), pts,
                               c.get("q_factor", 20.0))
    crit = min(env, key=lambda r: r["margin"])
    return {"case": c, "f_Hz": (omega / 2 / np.pi).tolist(), "D": lam["D"].tolist(), "critical": crit, "envelope": env}


def main():
    case = Path(sys.argv[1])
    r = run(case)
    c = r["critical"]
    print("f [Hz]:", " ".join(f"{f:.0f}" for f in r["f_Hz"]))
    print(f"min margin q_crit/q = {c['margin']:.2f} ({c['kind'] or 'stable'}, f={c['f_Hz']:.0f} Hz) at M={c['mach']:.3f}, "
          f"q={c['q']:.0f} Pa; constant-Mach velocity margin {np.sqrt(c['margin']):.2f}")
    sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    out = Path("results") / f"{case.stem}.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps({"commit": sha, **r}, indent=1, default=float))
    print("->", out)
