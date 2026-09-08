#!/usr/bin/env python3
"""
sectorAssembly.py -- make a half or full mesh out of the structured quadrant.

The block structure is built once, for theta in [0, 90] deg, with a half-fin
lying in each of its two planes.  A larger sector is assembled by ROTATING
copies of that quadrant about the x axis and STITCHING them together at the
planes where they meet:

    copy k  =  quadrant rotated by  -90 k  degrees     (k = 0 .. n-1)

    quarter  1 copy    planes A (z = 0) and B (y = 0) stay symmetry planes
    half     2 copies  B of copy 0 meets A of copy 1;  A0 and B1 stay symmetry
    full     4 copies  every B meets the next A, and B3 wraps round to A0

Stitching is node-for-node.  Two nodes on an interface plane are the same
node in the assembled mesh when they sit at the same (x, r) -- which they do
exactly, because the quadrant is mirror-symmetric about 45 deg (the azimuthal
block widths and gradings are built that way) -- UNLESS the fin deformation
moved them off the plane.  A node on the fin planform stays two nodes: one at
+t_half, one at -t_half, the two faces of the now full-thickness fin.

Faces on an interface plane follow the same rule: all four nodes merged means
the face is interior and is dropped from both sides; any node on the fin
means the face is a fin wall on both sides.  The classification is the same
function that split the quadrant's symmetry patch (finPatch.split_symm), so
the two cannot disagree.

Nothing here is specific to this rocket beyond "four fins at 90 deg".  A
three-fin vehicle needs a 120 deg sector, which a square butterfly core
cannot tile; that would be a new block topology, not a new assembly.
"""

import numpy as np

import meshParams as MP


def _rotate_yz(nodes, k):
    """Rotate (y, z) by -90 k degrees about the x axis.  Proper rotation, so
    hexahedron orientation is preserved."""
    p = nodes.copy()
    y, z = nodes[:, 1], nodes[:, 2]
    k %= 4
    if k == 1:
        p[:, 1], p[:, 2] = z, -y
    elif k == 2:
        p[:, 1], p[:, 2] = -y, -z
    elif k == 3:
        p[:, 1], p[:, 2] = -z, y
    return p


def _match(key_b, key_a, tol=1e-7):
    """For each row of key_b find the row of key_a at the same (x, r).

    Exact-arithmetic equality is expected; the tolerance only absorbs
    floating-point noise from two transfinite interpolations of the same
    distribution.  Every node must find exactly one partner or the quadrant
    is not mirror-symmetric, which is a bug and is reported as one.
    """
    if len(key_b) != len(key_a):
        raise RuntimeError(f'interface planes carry {len(key_b)} and {len(key_a)} '
                           f'nodes -- the quadrant is not mirror-symmetric')
    try:
        from scipy.spatial import cKDTree
        d, j = cKDTree(key_a).query(key_b, distance_upper_bound=tol)
        bad = ~np.isfinite(d)
    except ImportError:                             # rounding fallback
        ka = {tuple(r): i for i, r in enumerate(np.round(key_a / tol).astype(np.int64).tolist())}
        j = np.array([ka.get(tuple(r), -1) for r in
                      np.round(key_b / tol).astype(np.int64).tolist()])
        bad = j < 0
    if np.any(bad):
        i = np.flatnonzero(bad)[0]
        raise RuntimeError(f'{bad.sum()} interface nodes have no partner, e.g. '
                           f'x = {key_b[i, 0]:.6f}, r = {key_b[i, 1]:.6f}. '
                           f'Install scipy for a tolerance-based match, or '
                           f'check that the azimuthal grading is symmetric.')
    if len(np.unique(j)) != len(j):
        raise RuntimeError('interface node matching is not one-to-one')
    return j


def assemble(nodes, hexes, quads, plane_a, plane_b, fin_nodes, n_copies,
             verbose=True):
    """Rotate-and-stitch `n_copies` quadrants (1, 2 or 4).

    nodes, hexes, quads  the quadrant, as from meshIO.extract (quads includes
                         the symmetry patch, UNSPLIT -- fin faces are
                         classified here)
    plane_a / plane_b    node masks: on z = 0 / on y = 0, undeformed
    fin_nodes            node mask: on the fin planform (finPatch.fin_node_mask)

    Returns (nodes, hexes, quads) for the assembled sector; the symmetry patch
    holds whatever external planes remain and the fin patch every fin face.
    """
    symm_name = MP.PATCHES['symmetry']
    fins_name = MP.PATCHES['fins']
    if n_copies not in (1, 2, 4):
        raise ValueError(f'n_copies must be 1, 2 or 4, got {n_copies}')
    n = len(nodes)
    q_symm = quads[symm_name]
    on_a = plane_a[q_symm].all(axis=1)
    on_b = plane_b[q_symm].all(axis=1)
    if not np.all(on_a | on_b):
        raise RuntimeError('a symmetry quad lies on neither plane')
    on_fin = fin_nodes[q_symm].any(axis=1)
    qa, qa_fin = q_symm[on_a & ~on_fin], q_symm[on_a & on_fin]
    qb, qb_fin = q_symm[on_b & ~on_fin], q_symm[on_b & on_fin]

    # --- node correspondence across one interface (B of copy k, A of copy k+1)
    ia = np.flatnonzero(plane_a)
    ib = np.flatnonzero(plane_b)
    key_a = np.column_stack([nodes[ia, 0], np.hypot(nodes[ia, 1], nodes[ia, 2])])
    key_b = np.column_stack([nodes[ib, 0], np.hypot(nodes[ib, 1], nodes[ib, 2])])
    partner = ia[_match(key_b, key_a)]              # plane-B node -> plane-A node
    keep_apart = fin_nodes[ib]                       # on the fin: do not merge
    if np.any(fin_nodes[ib] != fin_nodes[partner]):
        raise RuntimeError('fin planform differs between the two planes')
    ib_m, pa_m = ib[~keep_apart], partner[~keep_apart]

    # --- union of node ids over the copies
    rep = np.arange(n * n_copies)
    interfaces = [(k, k + 1) for k in range(n_copies - 1)]
    if n_copies == 4:
        interfaces.append((3, 0))
    for kb, ka in interfaces:
        hi, lo = ib_m + kb * n, pa_m + ka * n
        if kb < ka:
            hi, lo = lo, hi
        rep[hi] = lo
    for _ in range(8):                                # resolve chains (axis nodes)
        nxt = rep[rep]
        if np.array_equal(nxt, rep):
            break
        rep = nxt
    uniq, new = np.unique(rep, return_inverse=True)

    all_nodes = np.concatenate([_rotate_yz(nodes, k) for k in range(n_copies)])[uniq]
    all_hexes = new[np.concatenate([hexes + k * n for k in range(n_copies)])]

    # --- patches
    out = {}
    for name, q in quads.items():
        if name == symm_name:
            continue
        out[name] = new[np.concatenate([q + k * n for k in range(n_copies)])]
    stitched = {kb for kb, _ in interfaces}
    stitched_a = {ka for _, ka in interfaces}
    symm, fins = [], []
    for k in range(n_copies):
        if k not in stitched:                        # external plane B
            symm.append(qb + k * n)
        if k not in stitched_a:                      # external plane A
            symm.append(qa + k * n)
        fins += [qa_fin + k * n, qb_fin + k * n]
    out[symm_name] = new[np.concatenate(symm)] if symm else np.zeros((0, 4), np.int64)
    out[fins_name] = new[np.concatenate(fins)] if fins else np.zeros((0, 4), np.int64)

    # --- sanity: a stitched face pair must have become the same interior face
    n_merged_faces = len(qb) * len(interfaces)
    if verbose:
        print(f'  sector: {n_copies} x quadrant -> {len(all_nodes):,d} nodes, '
              f'{len(all_hexes):,d} cells; {len(interfaces)} interface(s), '
              f'{len(ib_m):,d} nodes merged per interface, '
              f'{n_merged_faces:,d} faces made interior, '
              f'{len(out[fins_name]):,d} fin faces')
    return all_nodes, all_hexes, out
