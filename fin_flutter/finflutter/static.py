"""Static check: generalised forces of the CFD pressure on the modes, from postProcessing/surfaces/*/skin.vtp."""
from pathlib import Path

import numpy as np

from . import cfd


def shape(points, s, roll_deg, half_width=0.02):
    """w of every mode at points [m, 3]; the same interpolation as flutterModes.H."""
    fin, nodes, phi = s["fin"], s["nodes"], s["phi"]
    nx = len(np.unique(np.round(nodes[nodes[:, 1] == 0, 0], 12)))
    ny = len(nodes) // nx
    r = np.radians(roll_deg)
    er, en = np.array([0, np.cos(r), np.sin(r)]), np.array([0, -np.sin(r), np.cos(r)])
    p = np.asarray(points, float)
    rr = p * np.array([0, 1, 1])
    sp = rr @ er - fin.R
    inside = (sp >= -half_width) & (np.abs(rr @ en) <= half_width)
    eta = np.clip(sp / fin.span, 0, 1)
    c = fin.cr + eta * (fin.ct - fin.cr)
    xi = np.clip((p[:, 0] - fin.x_root_le - eta * fin.sweep) / c, 0, 1)
    fi, fj = xi * (nx - 1), eta * (ny - 1)
    i, j = np.minimum(fi.astype(int), nx - 2), np.minimum(fj.astype(int), ny - 2)
    a, b = (fi - i)[:, None], (fj - j)[:, None]
    W = phi.reshape(ny, nx, -1)
    w = (1 - a) * (1 - b) * W[j, i] + a * (1 - b) * W[j, i + 1] + (1 - a) * b * W[j + 1, i] + a * b * W[j + 1, i + 1]
    return w * inside[:, None], en


def skin(run):
    """(centres, area vectors, p) of the wall faces of the latest skin.vtp of a run."""
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy as v2n
    files = sorted((Path(run) / "postProcessing" / "surfaces").glob("*/skin.vtp"), key=lambda f: float(f.parent.name))
    if not files:
        raise FileNotFoundError(f"{run}: no postProcessing/surfaces/*/skin.vtp")
    rd = vtk.vtkXMLPolyDataReader(); rd.SetFileName(str(files[-1])); rd.Update()
    d = rd.GetOutput()
    nf = vtk.vtkPolyDataNormals(); nf.SetInputData(d); nf.ComputeCellNormalsOn(); nf.ComputePointNormalsOff()
    nf.SplittingOff(); nf.ConsistencyOff(); nf.AutoOrientNormalsOff(); nf.Update()
    cs = vtk.vtkCellSizeFilter(); cs.SetInputData(nf.GetOutput()); cs.ComputeAreaOn(); cs.Update()
    out = cs.GetOutput()
    cc = vtk.vtkCellCenters(); cc.SetInputData(out); cc.Update()
    X = v2n(cc.GetOutput().GetPoints().GetData())
    S = v2n(out.GetCellData().GetNormals()) * v2n(out.GetCellData().GetArray("Area"))[:, None]
    return X, S, v2n(out.GetCellData().GetArray("p"))


def generalised_forces(run, s, roll_deg):
    X, S, p = skin(run)
    flow = cfd.run_flow(run)
    w, en = shape(X, s, roll_deg)
    return -((p - flow["pInf"]) * (S @ en)) @ w, flow


def compare(run_alpha, run_zero, s):
    """CFD F_k(alpha) - F_k(0) against Ackeret piston theory for the same rigid incidence of fin 0."""
    roll = 45.0 if cfd.run_sector(run_alpha) == "half" else 0.0
    Fa, flow = generalised_forces(run_alpha, s, roll)
    F0, _ = generalised_forces(run_zero, s, roll)
    M, alpha = flow["Minf"], np.radians(flow["alpha"])
    q = 0.5 * flow["rhoInf"] * flow["Uinf"] ** 2
    theta = np.sin(alpha) * np.sin(np.radians(roll))
    ones = np.ones(len(s["nodes"]))
    integral = ones @ (s["mats"]["W"] @ s["phi"])
    piston = -q * 2 * (2 / np.sqrt(M * M - 1)) * theta * integral if M > 1 else np.full_like(integral, np.nan)
    return {"mach": M, "alpha_deg": np.degrees(alpha), "theta_deg": np.degrees(theta), "q": q,
            "F_cfd": Fa - F0, "F_piston": piston}
