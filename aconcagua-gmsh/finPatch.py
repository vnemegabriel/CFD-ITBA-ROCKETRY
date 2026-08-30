#!/usr/bin/env python3
"""
finPatch.py -- introduce the fins into an existing structured mesh.

Two steps, neither of which touches the block topology.

1. DEFORM the azimuthal coordinate:

       theta -> theta_f(x,r) + theta * (90 - 2 theta_f) / 90

   A node that was on the symmetry plane (theta = 0) lands at theta_f, i.e. at
   z = -r sin(theta_f) = -t_half(x,r), which IS the fin surface.  The geometry
   comes out exact rather than approximated, the mesh stays conforming, and not
   one block is added.  It works only because the real fin is BEVELLED: t goes
   to zero continuously at the leading and trailing edges, so the deformation
   relaxes to nothing there.  The tip is the one genuine discontinuity and is
   smeared over FIN_TIP_SMEAR.

   Both symmetry planes get it, because the quarter model contains two half
   fins -- one lying in each plane.

2. RETAG the faces that are now fin rather than symmetry.  A gmsh physical
   group is per-surface, and one block side face carries both, so the split has
   to happen at element level, in the .msh.
"""

import sys
import time

import numpy as np

import aconcaguaGeom as G
import meshParams as MP


def deform_nodes(gmsh, section=None, smear=None):
    """Apply the fin deformation to every node in the model."""
    section = section or MP.FIN_SECTION
    smear = MP.FIN_TIP_SMEAR if smear is None else smear
    moved = 0
    dmax = 0.0
    for dim, tag in gmsh.model.getEntities():
        t, c, par = gmsh.model.mesh.getNodes(dim, tag)
        if len(t) == 0:
            continue
        p = c.reshape(-1, 3).copy()
        x, y, z = p[:, 0], p[:, 1], p[:, 2]
        r = np.hypot(y, z)
        live = (r > 1e-9) & (x >= G.FIN_ROOT_LE - 1e-6) & (x <= G.FIN_TIP_TE + 1e-6)
        if not np.any(live):
            continue
        th = np.arctan2(-z[live], y[live])
        tf = G.fin_half_angle(x[live], r[live], section, smear)
        if not np.any(tf > 0):
            continue
        half = np.pi / 2.0
        thn = tf + th * (half - 2.0 * tf) / half
        rr = r[live]
        ny, nz = rr * np.cos(thn), -rr * np.sin(thn)
        d = np.hypot(ny - y[live], nz - z[live])
        dmax = max(dmax, float(d.max()))
        y[live], z[live] = ny, nz
        # gmsh 4.15 exposes only the singular setNode, so touch just the nodes
        # that actually moved -- ~1 % of the mesh, since the fin is local.
        idx = np.flatnonzero(live)[d > 1e-12]
        for k in idx:
            gmsh.model.mesh.setNode(int(t[k]), p[k].tolist(), [])
        moved += len(idx)
    return moved, dmax


# ------------------------------------------------------------ patch split ---
def split_symm(path, section=None, verbose=True):
    """Retag the symmetry-plane faces that are covered by the fin."""
    section = section or MP.FIN_SECTION
    t0 = time.time()
    src = open(path).read().split('\n')

    def block(name):
        i = src.index(f'${name}')
        j = src.index(f'$End{name}')
        return i, j

    # --- physical names ----------------------------------------------------
    i, j = block('PhysicalNames')
    names = src[i + 2:j]
    symm_tag = None
    used = set()
    for ln in names:
        f = ln.split(' ', 2)
        used.add(int(f[1]))
        if f[2].strip('"') == MP.PATCHES['symmetry']:
            symm_tag = int(f[1])
    if symm_tag is None:
        raise KeyError('no symmetry physical group in the mesh')
    fin_tag = max(used) + 1
    src[i + 1] = str(int(src[i + 1]) + 1)
    src[j:j] = [f'2 {fin_tag} "fins"']

    # --- nodes -------------------------------------------------------------
    i, j = block('Nodes')
    n = int(src[i + 1])
    arr = np.fromstring(' '.join(src[i + 2:j]), sep=' ').reshape(n, 4)
    lut = np.zeros(int(arr[:, 0].max()) + 1, np.int64)
    lut[arr[:, 0].astype(np.int64)] = np.arange(n)
    xyz = arr[:, 1:]

    # --- elements ----------------------------------------------------------
    i, j = block('Elements')
    out = []
    hit = 0
    pre = f' 3 2 {symm_tag} '
    for k in range(i + 2, j):
        ln = src[k]
        p = ln.find(pre)
        if p < 0:
            out.append(ln)
            continue
        f = ln.split()
        nodes = lut[np.array(f[5:9], dtype=np.int64)]
        c = xyz[nodes].mean(axis=0)
        r = float(np.hypot(c[1], c[2]))
        if G.fin_half_thickness(float(c[0]), r, section, tip_smear=0.0) > 1e-9:
            f[3] = str(fin_tag)
            hit += 1
            out.append(' '.join(f))
        else:
            out.append(ln)
    src[i + 2:j] = out

    open(path, 'w').write('\n'.join(src))
    if verbose:
        print(f'  fins patch: {hit:,d} faces retagged from '
              f'{MP.PATCHES["symmetry"]} in {time.time()-t0:.1f}s')
    return hit


if __name__ == '__main__':
    split_symm(sys.argv[1] if len(sys.argv) > 1 else 'aconcagua_body.msh')
