#!/usr/bin/env python3
"""
meshQuality.py -- OpenFOAM's own quality measures, computed here.

gmsh reports Jacobians and SICN; OpenFOAM's checkMesh rejects meshes on
non-orthogonality, skewness and face pyramids.  They are not the same numbers,
and a mesh that satisfies one can fail the other, so the honest thing is to
compute what the CONSUMER will compute -- before handing it over.

Definitions follow src/OpenFOAM/meshes/primitiveMesh/primitiveMeshCheck:
    non-orthogonality  angle between the face normal and the vector joining
                       the owner and neighbour cell centres
    skewness           distance from the face centre to where the owner->
                       neighbour line pierces the face, over the face size
"""

import numpy as np

# local node ordering of the 6 faces of a gmsh hexahedron (type 5)
HEX_FACES = np.array([[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4],
                      [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]])


def face_geometry_chunked(nodes, faces, chunk=1_500_000):
    """Face centres and area vectors, in chunks so a 14 M-face mesh fits in RAM."""
    N = len(faces)
    ctr = np.empty((N, 3)); area = np.empty((N, 3))
    for a in range(0, N, chunk):
        b = min(a + chunk, N)
        c, s = face_geometry(nodes[faces[a:b]])
        ctr[a:b] = c; area[a:b] = s
    return ctr, area


def face_geometry(pts):
    """Centre, area vector and 'size' of quad faces given as (n,4,3)."""
    c0 = pts.mean(axis=1)
    a = np.zeros((len(pts), 3))
    cw = np.zeros((len(pts), 3))
    tot = np.zeros(len(pts))
    for i in range(4):
        p, q = pts[:, i], pts[:, (i + 1) % 4]
        tri = np.cross(q - p, c0 - p) * 0.5
        mag = np.linalg.norm(tri, axis=1)
        a += tri
        cw += ((p + q + c0) / 3.0) * mag[:, None]
        tot += mag
    ok = tot > 1e-300
    ctr = np.where(ok[:, None], cw / np.where(ok, tot, 1.0)[:, None], c0)
    return ctr, a


def analyse(nodes, hexes, progress=None):
    """nodes: (N,3) float array.  hexes: (M,8) int index array."""
    M = len(hexes)
    faces = hexes[:, HEX_FACES.ravel()].reshape(M * 6, 4)
    owner = np.repeat(np.arange(M), 6)

    key = np.sort(faces, axis=1).astype(np.int32)
    view = np.ascontiguousarray(key).view([('a','<i4'),('b','<i4'),('c','<i4'),('d','<i4')]).ravel()
    order = np.argsort(view, kind='stable')
    ks = view[order]
    same = ks[1:] == ks[:-1]
    pair_i = order[:-1][same]
    pair_j = order[1:][same]
    del view, ks

    # every face must be used once (boundary) or twice (internal)
    used = np.zeros(M * 6, bool)
    used[pair_i] = used[pair_j] = True
    n_int, n_bnd = len(pair_i), int((~used).sum())
    del key
    triple = bool(np.any(np.isin(pair_i, pair_j)))

    # cell centres and volumes by decomposition about the cell centroid
    cc = nodes[hexes].mean(axis=1)
    fc, fa = face_geometry_chunked(nodes, faces)
    vol = np.zeros(M)
    np.add.at(vol, owner, np.einsum('ij,ij->i', fa, fc - cc[owner]) / 3.0)

    # orient every face outward of its owner, then measure the internal pairs
    sgn = np.sign(np.einsum('ij,ij->i', fa, fc - cc[owner]))
    fa = fa * sgn[:, None]

    magS_all = np.linalg.norm(fa, axis=1)
    o, n = owner[pair_i], owner[pair_j]
    Sf = fa[pair_i]
    Cf = fc[pair_i]
    d = cc[n] - cc[o]
    magS = np.linalg.norm(Sf, axis=1)
    magd = np.linalg.norm(d, axis=1)
    cosa = np.clip(np.einsum('ij,ij->i', Sf, d) / (magS * magd), -1.0, 1.0)
    nonortho = np.degrees(np.arccos(cosa))

    # skewness: where o->n pierces the face plane, versus the face centre
    nrm = Sf / magS[:, None]
    t = np.einsum('ij,ij->i', Cf - cc[o], nrm) / np.einsum('ij,ij->i', d, nrm)
    pierce = cc[o] + t[:, None] * d
    skew = np.linalg.norm(pierce - Cf, axis=1) / np.sqrt(magS)

    # aspect ratio: OpenFOAM's definition, 1/6 * sum|Sf| * cell size / volume
    sumS = np.zeros(M)
    np.add.at(sumS, owner, magS_all)
    with np.errstate(divide='ignore', invalid='ignore'):
        ar = (1.0 / 6.0) * sumS * np.cbrt(np.abs(vol)) / np.abs(vol)

    worst = np.argsort(nonortho)[-40:][::-1]
    hot = np.column_stack([Cf[worst], nonortho[worst]])
    bad = nonortho > 40.0
    return dict(hot=hot, hot_all=np.column_stack([Cf[bad], nonortho[bad]]),
                n_cells=M, n_internal=n_int, n_boundary=n_bnd,
                triple_face=bool(triple),
                vol_min=float(vol.min()), vol_max=float(vol.max()),
                vol_total=float(vol.sum()), n_negative=int((vol <= 0).sum()),
                nonortho_max=float(nonortho.max()),
                nonortho_mean=float(nonortho.mean()),
                nonortho_gt70=int((nonortho > 70).sum()),
                nonortho_gt40=int((nonortho > 40).sum()),
                skew_max=float(skew.max()), skew_mean=float(skew.mean()),
                skew_gt4=int((skew > 4).sum()),
                ar_max=float(np.nanmax(ar)), ar_mean=float(np.nanmean(ar)))
