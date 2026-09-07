#!/usr/bin/env python3
"""
meshIO.py -- get the mesh OUT of gmsh as arrays, and write it for OpenFOAM.

After `gmsh.model.mesh.generate(3)` everything downstream (wall snap, fin
deformation, sector assembly, quality audit, export) works on plain numpy
arrays:

    nodes   (n, 3) float   x, y, z
    hexes   (m, 8) int     0-based node indices, gmsh hexahedron ordering
    quads   {patch: (k, 4) int}   boundary faces, by patch name

Doing it this way rather than through gmsh's per-node `setNode` API is both
faster (one vectorised pass instead of a Python loop over every moved node)
and simpler: the .msh writer here emits exactly the v2.2 ASCII layout that
gmshToFoam consumes, with the physical names it will turn into patches, and
nothing else.  Stray spline control points never reach the file because only
nodes referenced by a hexahedron are written.
"""

import numpy as np

# gmsh element type codes
QUAD4, HEX8 = 3, 5


# ----------------------------------------------------------------- extract --
def extract(bodymesh):
    """Pull nodes / hexes / boundary quads from the live gmsh model.

    `bodymesh` is a built BodyMesh (buildHexBody.py); its Deck knows which
    surfaces carry which patch, and which surfaces are the body wall.
    Returns a dict with compact 0-based arrays plus the index sets of the
    lateral-wall and nose-cap nodes (for the snap step).
    """
    import gmsh
    nt, nc, _ = gmsh.model.mesh.getNodes()
    nodes = nc.reshape(-1, 3)
    lut = np.full(int(nt.max()) + 1, -1, np.int64)
    lut[nt.astype(np.int64)] = np.arange(len(nt))

    _, en = gmsh.model.mesh.getElementsByType(HEX8)
    hexes = lut[en.astype(np.int64)].reshape(-1, 8)

    quads, owner = {}, {}
    for name, tags in bodymesh.D.surf_tags.items():
        parts = []
        for s in sorted(set(tags)):
            if s in owner and owner[s] != name:
                raise RuntimeError(f'surface {s} is in patches {owner[s]} and {name}')
            owner[s] = name
            _, qn = gmsh.model.mesh.getElementsByType(QUAD4, s)
            if len(qn):
                parts.append(lut[qn.astype(np.int64)].reshape(-1, 4))
        quads[name] = np.concatenate(parts) if parts else np.zeros((0, 4), np.int64)

    def surf_nodes(surfs):
        out = []
        for s in sorted(surfs):
            t, _, _ = gmsh.model.mesh.getNodes(2, s, includeBoundary=True)
            out.append(lut[t.astype(np.int64)])
        return np.unique(np.concatenate(out)) if out else np.zeros(0, np.int64)

    lateral = surf_nodes(bodymesh.wall_lateral)
    cap = surf_nodes(bodymesh.wall_cap)

    # keep only nodes a hexahedron references: drops spline control points,
    # which gmsh would otherwise export as unreferenced nodes
    used = np.unique(hexes)
    remap = np.full(len(nodes), -1, np.int64)
    remap[used] = np.arange(len(used))
    nodes = nodes[used]
    hexes = remap[hexes]
    for k in quads:
        q = remap[quads[k]]
        if np.any(q < 0):
            raise RuntimeError(f'patch {k} references a node no cell uses')
        quads[k] = q
    lateral = remap[lateral]; lateral = lateral[lateral >= 0]
    cap = remap[cap]; cap = cap[cap >= 0]
    return dict(nodes=nodes, hexes=hexes, quads=quads,
                lateral_idx=lateral, cap_idx=cap, n_dropped=int(len(remap) - len(used)))


# ------------------------------------------------------------------- write --
def _write_rows(f, arr, fmt, chunk=400_000):
    """Write integer/float rows with a fixed printf format, in chunks."""
    for a in range(0, len(arr), chunk):
        rows = arr[a:a + chunk].tolist()
        f.write('\n'.join(fmt % tuple(r) for r in rows))
        f.write('\n')


def write_msh2(path, nodes, hexes, quads, patch_order=None, volume_name='internal'):
    """Write a gmsh v2.2 ASCII file that gmshToFoam turns into patches.

    Patches are numbered in `patch_order` (default: dict order); the volume
    physical group comes last.  Element ids are 1-based and contiguous, quads
    first.  Returns the {patch: physical tag} map.
    """
    names = list(patch_order or quads.keys())
    names = [n for n in names if len(quads[n])]
    tags = {n: i + 1 for i, n in enumerate(names)}
    vtag = len(names) + 1

    with open(path, 'w', newline='\n') as f:
        f.write('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n')
        f.write(f'$PhysicalNames\n{len(names) + 1}\n')
        for n in names:
            f.write(f'2 {tags[n]} "{n}"\n')
        f.write(f'3 {vtag} "{volume_name}"\n$EndPhysicalNames\n')

        f.write(f'$Nodes\n{len(nodes)}\n')
        ids = np.arange(1, len(nodes) + 1)
        _write_rows(f, np.column_stack([ids, nodes]), '%d %.16g %.16g %.16g')
        f.write('$EndNodes\n')

        n_el = sum(len(quads[n]) for n in names) + len(hexes)
        f.write(f'$Elements\n{n_el}\n')
        eid = 1
        for n in names:
            q = quads[n]
            k = len(q)
            block = np.column_stack([np.arange(eid, eid + k), np.full(k, QUAD4),
                                     np.full(k, 2), np.full(k, tags[n]),
                                     np.full(k, tags[n]), q + 1])
            _write_rows(f, block, '%d %d %d %d %d %d %d %d %d')
            eid += k
        k = len(hexes)
        block = np.column_stack([np.arange(eid, eid + k), np.full(k, HEX8),
                                 np.full(k, 2), np.full(k, vtag), np.full(k, vtag),
                                 hexes + 1])
        _write_rows(f, block, '%d %d %d %d %d %d %d %d %d %d %d %d %d')
        f.write('$EndElements\n')
    return tags


# ---------------------------------------------------------------- sidecars --
def _foam(v):
    if isinstance(v, bool):
        return 'true' if v else 'false'
    if isinstance(v, (list, tuple)):
        return '(' + ' '.join(_foam(x) for x in v) + ')'
    if isinstance(v, float):
        return f'{v:.8g}'
    return str(v)


def write_meshinfo(path, info):
    """An OpenFOAM dictionary describing the mesh, for the case to #include.

    The case reads Aref / lRef / the wall patch list / the patch types from
    here, so the mesh is the single source of truth for what the case must
    agree with.  `info` keys map straight to dictionary entries; the
    `patchTypes` entry is a sub-dictionary.
    """
    lines = ['/*--------------------------------*- C++ -*----------------------------------*\\',
             '|  Written by mesh/build.py -- describes the .msh this case was built from.   |',
             '|  The case #includes it: do not edit by hand, rebuild the mesh instead.      |',
             '\\*---------------------------------------------------------------------------*/',
             'FoamFile', '{', '    version     2.0;', '    format      ascii;',
             '    class       dictionary;', '    object      meshInfo;', '}', '']
    for k, v in info.items():
        if k == 'patchTypes':
            lines.append('patchTypes')
            lines.append('{')
            for p, t in v.items():
                lines.append(f'    {p:<16}{t};')
            lines.append('}')
        else:
            lines.append(f'{k:<20}{_foam(v)};')
    lines.append('')
    open(path, 'w', newline='\n').write('\n'.join(lines))


def write_params(path, params, header=''):
    """The full parameter set as a Python preset -- feed it back to build.py
    with --preset to rebuild exactly this mesh."""
    out = ['# Parameter record written by build.py.  Reusable as a preset:',
           '#     python build.py --preset <this file>', '']
    if header:
        out += ['# ' + l for l in header.splitlines()] + ['']
    for k, v in params.items():
        out.append(f'{k} = {v!r}')
    out.append('')
    open(path, 'w', newline='\n').write('\n'.join(out))
