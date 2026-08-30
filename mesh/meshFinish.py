#!/usr/bin/env python3
"""
meshFinish.py -- build, snap the wall onto the analytic surface, audit, export.

The snap step is the one place where this workflow does anything resembling
projection -- and the point of measuring the deviation FIRST is to show how
little there is to do.  gmsh's transfinite surfaces interpolate between exact
boundary curves, so the interior of a wall patch is already within a fraction
of a micron of the true surface of revolution.  Compare with a snap phase that
has to move a background hex corner up to a whole cell to reach the STL, and
can give up part way with `Displacement scaling for error reduction set to 0`.
"""

import os
import sys
import time

import numpy as np
import gmsh

import aconcaguaGeom as G
import meshParams as MP
from buildHexBody import BodyMesh
from meshQuality import analyse
import finPatch


def wall_deviation(node_xyz):
    """Radial distance of each node from the exact surface of revolution."""
    x = node_xyz[:, 0]
    r = np.hypot(node_xyz[:, 1], node_xyz[:, 2])
    return r - np.asarray(G.r_body(np.clip(x, 0.0, G.X_BASE)))


def snap(mesh, verbose=True):
    """Move every wall node exactly onto the analytic surface."""
    report = {}

    # --- lateral walls: correct the RADIUS, keep x and theta ---------------
    tags, coord = [], []
    for s in sorted(mesh.wall_lateral):
        t, c, _ = gmsh.model.mesh.getNodes(2, s, includeBoundary=True)
        tags.append(t); coord.append(c.reshape(-1, 3))
    if tags:
        t = np.concatenate(tags); c = np.concatenate(coord)
        t, idx = np.unique(t, return_index=True); c = c[idx]
        dev = wall_deviation(c)
        report['lateral_before'] = float(np.abs(dev).max())
        rr = np.hypot(c[:, 1], c[:, 2])
        scale = np.where(rr > 1e-14, (rr - dev) / np.where(rr > 1e-14, rr, 1.0), 1.0)
        c[:, 1] *= scale; c[:, 2] *= scale
        for tag, p in zip(t, c):
            gmsh.model.mesh.setNode(int(tag), p.tolist(), [])
        report['lateral_after'] = float(np.abs(wall_deviation(c)).max())
        report['lateral_nodes'] = int(len(t))

    # --- the butterfly cap: correct X, keep (y, z) -------------------------
    tags, coord = [], []
    for s in sorted(mesh.wall_cap):
        t, c, _ = gmsh.model.mesh.getNodes(2, s, includeBoundary=True)
        tags.append(t); coord.append(c.reshape(-1, 3))
    if tags:
        t = np.concatenate(tags); c = np.concatenate(coord)
        t, idx = np.unique(t, return_index=True); c = c[idx]
        report['cap_before'] = float(np.abs(wall_deviation(c)).max())
        rr = np.hypot(c[:, 1], c[:, 2])
        c[:, 0] = np.where(rr > 1e-12, G.x_nose_of_r(np.maximum(rr, 1e-12)), 0.0)
        for tag, p in zip(t, c):
            gmsh.model.mesh.setNode(int(tag), p.tolist(), [])
        report['cap_after'] = float(np.abs(wall_deviation(c)).max())
        report['cap_nodes'] = int(len(t))

    if verbose:
        print('\n  wall projection (max |r - r_exact|)')
        for k in ('lateral', 'cap'):
            if f'{k}_before' in report:
                print(f"    {k:<9} {report[f'{k}_nodes']:>8,d} nodes   "
                      f"before {report[f'{k}_before']*1e9:9.2f} nm   "
                      f"after {report[f'{k}_after']*1e9:8.2f} nm")
    return report


def hex_arrays():
    """Node coordinates and hex connectivity as dense 0-based numpy arrays."""
    nt, nc, _ = gmsh.model.mesh.getNodes()
    nc = nc.reshape(-1, 3)
    lut = np.zeros(int(nt.max()) + 1, np.int32)
    lut[nt.astype(np.int64)] = np.arange(len(nt), dtype=np.int32)
    et, en = gmsh.model.mesh.getElementsByType(5)
    hx = lut[en.astype(np.int64)].reshape(-1, 8)
    return nc, hx


def main():
    t0 = time.time()
    coarse = os.environ.get('AG_COARSE', '1')
    print(f'building  (AG_COARSE={coarse})')
    m = BodyMesh()
    m.build()
    print('  geometry:', m.D.stats())
    gmsh.model.mesh.generate(3)
    print(f'  meshed in {time.time()-t0:.1f}s')

    n_stray = m.D.drop_control_nodes()
    rep = snap(m)

    if MP.FINS_ON:
        nmv, dmx = finPatch.deform_nodes(gmsh)
        print(f'\n  fins: {nmv:,d} nodes deformed, max displacement '
              f'{dmx*1e3:.3f} mm  (section {MP.FIN_SECTION!r})')

    nodes, hexes = hex_arrays()
    print(f'\n  nodes {len(nodes):,d}   hexes {len(hexes):,d}   (dropped {n_stray:,d} spline control nodes)')

    q = analyse(nodes, hexes)
    print(f"""
  ---- mesh audit -------------------------------------------------
  cells                       {q['n_cells']:>14,d}   (100% hexahedra)
  internal faces              {q['n_internal']:>14,d}
  boundary faces              {q['n_boundary']:>14,d}
  a face shared by 3+ cells   {str(q['triple_face']):>14}
  cells with volume <= 0      {q['n_negative']:>14,d}
  total volume                {q['vol_total']:>14,.3f} m3
  min / max cell volume       {q['vol_min']:>10.3e} / {q['vol_max']:.3e}
  non-orthogonality  max      {q['nonortho_max']:>14.2f} deg
                     mean     {q['nonortho_mean']:>14.2f} deg
                     > 70 deg {q['nonortho_gt70']:>14,d} faces
                     > 40 deg {q['nonortho_gt40']:>14,d} faces
  skewness           max      {q['skew_max']:>14.3f}
                     mean     {q['skew_mean']:>14.3f}
  aspect ratio       max      {q['ar_max']:>14.1f}
                     mean     {q['ar_mean']:>14.2f}
  -----------------------------------------------------------------""")

    if os.environ.get('AG_WRITE', '1') == '1':
        gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
        gmsh.option.setNumber('Mesh.Binary', 0)
        gmsh.write('aconcagua_body.msh')
        print('  wrote aconcagua_body.msh  (v2.2 ASCII, for gmshToFoam)')
        if MP.FINS_ON:
            finPatch.split_symm('aconcagua_body.msh')
    gmsh.finalize()
    print(f'  total {time.time()-t0:.1f}s')
    return q


if __name__ == '__main__':
    main()
