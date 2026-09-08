#!/usr/bin/env python3
"""
meshFinish.py -- the mesh pipeline, from parameters to arrays.

    generate_quadrant()   build blocks -> gmsh mesh -> arrays -> snap wall
                          -> fin edge fit -> fin deformation -> fin faces
                          classified
    audit()               OpenFOAM's quality measures on any (nodes, hexes)

build.py is the command-line front end; this module is what it calls.  Use
it directly from Python when you want the arrays (e.g. to plot a cut, or to
compare two parameter sets in one session).

The snap step is the one place where this workflow does anything resembling
projection -- and the point of measuring the deviation FIRST is to show how
little there is to do.  gmsh's transfinite surfaces interpolate between exact
boundary curves, so the interior of a wall patch is already within a fraction
of a micron of the true surface of revolution.  Compare with a snap phase that
has to move a background hex corner up to a whole cell to reach an STL, and
can give up part way.
"""

import time

import numpy as np

import aconcaguaGeom as G
import meshParams as MP
import finPatch
import finEdge
import meshIO
from meshQuality import analyse

PLANE_TOL = 1e-12


def wall_deviation(xyz):
    """Radial distance of each node from the exact surface of revolution."""
    x = xyz[:, 0]
    r = np.hypot(xyz[:, 1], xyz[:, 2])
    return r - np.asarray(G.r_body(np.clip(x, 0.0, G.X_BASE)))


def snap(nodes, lateral_idx, cap_idx, verbose=True):
    """Move every wall node exactly onto the analytic surface, in place."""
    rep = {}
    if len(lateral_idx):                    # lateral walls: fix r, keep x and theta
        c = nodes[lateral_idx]
        dev = wall_deviation(c)
        rep['lateral_before'] = float(np.abs(dev).max())
        rr = np.hypot(c[:, 1], c[:, 2])
        scale = np.where(rr > 1e-14, (rr - dev) / np.where(rr > 1e-14, rr, 1.0), 1.0)
        c[:, 1] *= scale
        c[:, 2] *= scale
        nodes[lateral_idx] = c
        rep['lateral_after'] = float(np.abs(wall_deviation(c)).max())
        rep['lateral_nodes'] = int(len(lateral_idx))
    if len(cap_idx):                        # nose cap: fix x, keep (y, z)
        c = nodes[cap_idx]
        rep['cap_before'] = float(np.abs(wall_deviation(c)).max())
        rr = np.hypot(c[:, 1], c[:, 2])
        c[:, 0] = np.where(rr > 1e-12, G.x_nose_of_r(np.maximum(rr, 1e-12)), 0.0)
        nodes[cap_idx] = c
        rep['cap_after'] = float(np.abs(wall_deviation(c)).max())
        rep['cap_nodes'] = int(len(cap_idx))
    if verbose:
        print('\n  wall projection (max |r - r_exact|)')
        for k in ('lateral', 'cap'):
            if f'{k}_before' in rep:
                print(f"    {k:<9} {rep[f'{k}_nodes']:>8,d} nodes   "
                      f"before {rep[f'{k}_before']*1e9:9.2f} nm   "
                      f"after {rep[f'{k}_after']*1e9:8.2f} nm")
    return rep


def generate_quadrant(verbose=True):
    """Parameters -> the meshed, snapped, fin-deformed quadrant as arrays.

    Returns a dict:
        nodes, hexes        the quadrant
        quads               {patch: faces}; the symmetry patch is UNSPLIT
        plane_a, plane_b    node masks on z = 0 and y = 0 (undeformed)
        fin_nodes           node mask on the fin planform
        stats, snap, fins   reports
    """
    import gmsh
    from buildHexBody import BodyMesh
    t0 = time.time()
    m = BodyMesh()
    m.build()
    stats = m.D.stats()
    if verbose:
        print('  geometry:', stats, ' blocks:', len(m.segs))
    gmsh.model.mesh.generate(3)
    if verbose:
        print(f'  meshed in {time.time()-t0:.1f}s')
    q = meshIO.extract(m)
    gmsh.finalize()

    nodes = q['nodes']
    rep_snap = snap(nodes, q['lateral_idx'], q['cap_idx'], verbose)

    # Bend the axial stations onto the fin leading and trailing edges, so the
    # face classification below has an exact outline to work from instead of a
    # one-cell staircase.  AFTER snap (which is what pins the wall the shear
    # then leaves alone) and BEFORE the fin deformation, which reads (x, r).
    edge = dict(on=False, moved=0, dmax=0.0)
    if MP.FINS_ON and MP.fin_edge_fit():
        rep_edge = finEdge.check()
        nodes, edge['moved'], edge['dmax'] = finEdge.warp(nodes)
        edge.update(on=True, min_stretch=rep_edge['min_stretch'])
        if verbose:
            print(f"\n  fin edges: {edge['moved']:,d} nodes sheared, max shift "
                  f"{edge['dmax']*1e3:.1f} mm, worst block stretch "
                  f"{rep_edge['min_stretch']:.2f}")

    # plane membership BEFORE the deformation moves fin nodes off the planes
    plane_a = np.abs(nodes[:, 2]) < PLANE_TOL
    plane_b = np.abs(nodes[:, 1]) < PLANE_TOL

    fins = dict(on=bool(MP.FINS_ON), moved=0, dmax=0.0)
    fin_nodes = np.zeros(len(nodes), bool)
    if MP.FINS_ON:
        fin_nodes = finPatch.fin_node_mask(nodes)
        nodes, fins['moved'], fins['dmax'] = finPatch.deform(nodes)
        if verbose:
            print(f"\n  fins: {fins['moved']:,d} nodes deformed, max displacement "
                  f"{fins['dmax']*1e3:.3f} mm  (section {MP.FIN_SECTION!r})")

    if verbose:
        print(f"\n  quadrant: {len(nodes):,d} nodes   {len(q['hexes']):,d} hexes   "
              f"(dropped {q['n_dropped']:,d} spline control nodes)")
    return dict(nodes=nodes, hexes=q['hexes'], quads=q['quads'],
                plane_a=plane_a, plane_b=plane_b, fin_nodes=fin_nodes,
                stats=stats, snap=rep_snap, fins=fins, fin_edge=edge,
                seconds=time.time() - t0)


def audit(nodes, hexes, quads, verbose=True, label='mesh audit'):
    """Run OpenFOAM's checks here, before OpenFOAM does.

    Also verifies that every boundary face of the cell set is in exactly one
    patch -- the check that nothing lands in defaultFaces and, after a sector
    assembly, that every stitched face really did become interior.
    """
    q = analyse(nodes, hexes)
    q['n_patch_faces'] = int(sum(len(v) for v in quads.values()))
    q['patches_consistent'] = q['n_patch_faces'] == q['n_boundary']
    if verbose:
        print(f"""
  ---- {label} {'-' * (60 - len(label))}
  cells                       {q['n_cells']:>14,d}   (100% hexahedra)
  internal faces              {q['n_internal']:>14,d}
  boundary faces              {q['n_boundary']:>14,d}""")
        print(f"  faces in patches            {q['n_patch_faces']:>14,d}   "
              f"{'OK' if q['patches_consistent'] else '!! MISMATCH -- defaultFaces'}")
        print(f"""  a face shared by 3+ cells   {str(q['triple_face']):>14}
  cells with volume <= 0      {q['n_negative']:>14,d}
  total volume                {q['vol_total']:>14,.3f} m3
  min / max cell volume       {q['vol_min']:>10.3e} / {q['vol_max']:.3e}
  non-orthogonality  max      {q['nonortho_max']:>14.2f} deg
                     mean     {q['nonortho_mean']:>14.2f} deg
                     > 70 deg {q['nonortho_gt70']:>14,d} faces
                     > 40 deg {q['nonortho_gt40']:>14,d} faces
  skewness           max      {q['skew_max']:>14.3f}   (internal {q['skew_max_internal']:.2f}, boundary {q['skew_max_boundary']:.2f}; > 4: {q['skew_gt4']:,d})
                     mean     {q['skew_mean']:>14.3f}
  aspect ratio       max      {q['ar_max']:>14.1f}
                     mean     {q['ar_mean']:>14.2f}
  -----------------------------------------------------------------""")
    return q
