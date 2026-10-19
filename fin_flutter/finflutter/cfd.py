"""Exchange with the OpenFOAM cases of flightModel: constant/flutterModes, run.info, genForces.dat."""
import math
import re
from pathlib import Path

import numpy as np


def pulse(nmodes, mode, amplitude, maxabs, fmax):
    """Gaussian pulse of one mode: eta(t) = coef exp(-((t - t0)/tau)^2), spectrum 1 % of its peak at fmax."""
    tau = 2 * math.sqrt(math.log(100)) / (2 * math.pi * fmax)
    coef = np.zeros(nmodes)
    coef[mode - 1] = amplitude / maxabs[mode - 1]
    return {"type": "pulse", "amplitude": amplitude, "coef": coef, "t0": 3.5 * tau, "tau": tau}


def harmonic(phi, v, amplitude, omega, n_ramp=2):
    """eta_k(t) = coef_k ramp(t) sin(omega t + phase_k), so that max |sum v_k phi_k| = amplitude."""
    v = np.asarray(v, complex)
    peak = np.abs(phi @ v).max()
    return {"type": "harmonic", "amplitude": amplitude, "coef": amplitude * np.abs(v) / peak, "phase": np.angle(v),
            "omega": omega, "tRamp": n_ramp * 2 * math.pi / omega}


def no_signal(nmodes):
    return {"type": "none", "coef": np.zeros(nmodes)}


def _list(a):
    return "(" + " ".join(f"{x:.10g}" for x in a) + ")"


def write_flutter_modes(path, s, signal, roll_deg=0.0, half_width=0.02):
    """OpenFOAM dictionary read by flightModel/case-euler/system/flutterModes.H."""
    fin, nodes, phi, omega = s["fin"], s["nodes"], s["phi"], s["omega"]
    nx = len(np.unique(np.round(nodes[nodes[:, 1] == 0, 0], 12)))
    ny = len(nodes) // nx
    sig = "\n".join(f"    {k:<10}{_list(v) if np.ndim(v) else f'{v:.10g}' if not isinstance(v, str) else v};"
                    for k, v in signal.items())
    text = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      flutterModes;
}}
// Written by fin_flutter (python -m finflutter modes).  Modes of fin 0, the fin
// at rollDeg about +x; w [m] per unit modal coordinate, mass-normalised, on the
// FE nodes: index j*nXi + i, xi = i/(nXi-1) along the chord, eta = j/(nEta-1)
// along the span.

xRootLE     {fin.x_root_le:.10g};
rBody       {fin.R:.10g};
span        {fin.span:.10g};
rootChord   {fin.cr:.10g};
tipChord    {fin.ct:.10g};
sweep       {fin.sweep:.10g};
rollDeg     {roll_deg:.10g};
halfWidth   {half_width:.10g};
nXi         {nx};
nEta        {ny};
nModes      {phi.shape[1]};
omega       {_list(omega)};
maxAbs      {_list(np.abs(phi).max(axis=0))};

modes
(
{chr(10).join(_list(phi[:, k]) for k in range(phi.shape[1]))}
);

signal
{{
{sig}
}}
"""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text)
    return Path(path)


def read_foam_scalars(path):
    """Plain 'name value;' entries of an OpenFOAM dictionary (no #include, no #eval)."""
    out = {}
    for m in re.finditer(r"^\s*(\w+)\s+([-+0-9.eE]+)\s*;", Path(path).read_text(), re.M):
        out[m.group(1)] = float(m.group(2))
    return out


def run_flow(run):
    """Free-stream values written by newCase.sh in run.info (block flow)."""
    text = (Path(run) / "run.info").read_text()
    block = text[text.index("\nflow"):]
    block = block[:block.index("\n}")]
    return {m.group(1): float(m.group(2)) for m in re.finditer(r"^\s*(\w+)\s+([-+0-9.eE]+)\s*;", block, re.M)}


def run_sector(run):
    info = Path(run) / "constant" / "meshInfo"
    m = re.search(r"^\s*sector\s+(\w+)\s*;", info.read_text(), re.M)
    return m.group(1) if m else "full"


def read_genforces(run):
    """(t, eta[nt, n], F[nt, n]) from postProcessing/genForces/*/genForces.dat, all start times joined."""
    files = sorted((Path(run) / "postProcessing" / "genForces").glob("*/genForces.dat"), key=lambda f: float(f.parent.name))
    if not files:
        raise FileNotFoundError(f"{run}: no postProcessing/genForces/*/genForces.dat")
    d = np.vstack([np.loadtxt(f, comments="#", ndmin=2) for f in files])
    _, keep = np.unique(d[:, 0], return_index=True)
    d = d[keep]
    n = (d.shape[1] - 1) // 2
    return d[:, 0], d[:, 1:1 + n], d[:, 1 + n:]
