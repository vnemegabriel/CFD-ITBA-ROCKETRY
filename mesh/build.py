#!/usr/bin/env python3
"""
build.py -- the one command that makes a mesh.

    python build.py --plan                       print the cell plan, build nothing
    python build.py                              defaults: quarter, H_SCALE 1  (~3.5 M cells)
    python build.py --preset coarse              presets/coarse.py on top of the defaults
    python build.py --scale 2 --sector half      2x coarser, 180 deg with one symmetry plane
    python build.py --sector full --no-fins      360 deg, clean body of revolution
    python build.py --set ZONE_R=[0.28,0.35,2,11.775] --set ZONE_H=[0.006,0.02,0.2,1.6] --set WAKE_ZONE_K=1
    python build.py --preset output/aconcagua_half_s2.params.py     rebuild exactly that mesh

Precedence, lowest to highest: meshParams.py defaults < --preset < --scale /
--sector / --fins < --set.  The final parameter set is written next to the
.msh as <name>.params.py, and an OpenFOAM dictionary <name>.meshInfo tells
the case what it is getting (sector, Aref, patch types).  case-*/Allmesh reads
both.
"""

import argparse
import ast
import os
import runpy
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import aconcaguaGeom as G          # noqa: E402
import meshParams as MP            # noqa: E402


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--plan', action='store_true', help='print the cell plan and exit')
    ap.add_argument('--preset', action='append', default=[], metavar='NAME|FILE',
                    help='a file in presets/ (by name) or any .py of NAME = value lines; repeatable')
    ap.add_argument('--scale', type=float, metavar='H', help='H_SCALE: multiplies every cell size')
    ap.add_argument('--sector', choices=sorted(MP.SECTORS), help='quarter (default), half or full')
    g = ap.add_mutually_exclusive_group()
    g.add_argument('--fins', dest='fins', action='store_true', default=None)
    g.add_argument('--no-fins', dest='fins', action='store_false')
    ap.add_argument('--set', action='append', default=[], metavar='KEY=VALUE',
                    help='override any parameter; VALUE is a Python literal')
    ap.add_argument('--out', metavar='FILE.msh',
                    help='output path (default output/<auto name>.msh).  Worth pointing '
                         'outside a synced folder: OneDrive uploading a 2 GB .msh as it '
                         'is written can double the time')
    ap.add_argument('--name', metavar='NAME', help='basename for the output files (no extension)')
    ap.add_argument('--no-write', action='store_true', help='build and audit, write nothing')
    ap.add_argument('--no-audit', action='store_true',
                    help='skip the quality audit of the ASSEMBLED mesh (the quadrant is always '
                         'audited).  The audit of a 14 M-cell full mesh needs about 10 GB')
    ap.add_argument('--quiet', action='store_true')
    return ap.parse_args(argv)


def load_preset(spec):
    """A preset is a Python file of NAME = value assignments."""
    path = spec
    if not os.path.exists(path):
        cand = os.path.join(HERE, 'presets', spec if spec.endswith('.py') else spec + '.py')
        if not os.path.exists(cand):
            have = sorted(f[:-3] for f in os.listdir(os.path.join(HERE, 'presets')) if f.endswith('.py'))
            raise SystemExit(f'preset {spec!r} not found. Available: {have}, or give a file path.')
        path = cand
    ns = runpy.run_path(path)
    values = {k: v for k, v in ns.items() if k.isupper() and not k.startswith('_')}
    MP.apply_overrides(values, source=path)
    return path, values


def apply_cli(a):
    """Layer the command-line choices onto meshParams.  Returns a provenance list."""
    prov = []
    for p in a.preset:
        path, vals = load_preset(p)
        prov.append(f'preset {path}: {", ".join(f"{k}={v!r}" for k, v in vals.items())}')
    direct = {}
    if a.scale is not None:
        direct['H_SCALE'] = a.scale
    if a.sector is not None:
        direct['SECTOR'] = a.sector
    if a.fins is not None:
        direct['FINS_ON'] = a.fins
    if direct:
        MP.apply_overrides(direct, source='command line')
        prov.append('flags: ' + ', '.join(f'{k}={v!r}' for k, v in direct.items()))
    sets = {}
    for kv in a.set:
        if '=' not in kv:
            raise SystemExit(f'--set expects KEY=VALUE, got {kv!r}')
        k, v = kv.split('=', 1)
        try:
            sets[k.strip()] = ast.literal_eval(v.strip())
        except (ValueError, SyntaxError):
            sets[k.strip()] = v.strip()          # a bare string such as wedge
    if sets:
        MP.apply_overrides(sets, source='--set')
        prov.append('--set: ' + ', '.join(f'{k}={v!r}' for k, v in sets.items()))
    return prov


def auto_name():
    s = f'{MP.H_SCALE:g}'.replace('.', 'p')
    return f'aconcagua_{MP.SECTOR}_s{s}' + ('' if MP.FINS_ON else '_nofins')


def mesh_info(n_cells, d, quads):
    """What the OpenFOAM case needs to know about this mesh."""
    P = MP.PATCHES
    frac = MP.symmetry_fraction()
    ref = G.reference_values(frac)
    types = {P['inlet']: 'patch', P['outlet']: 'patch', P['farfield']: 'patch',
             P['symmetry']: 'symmetry', P['nose']: 'wall', P['body']: 'wall',
             P['tail']: 'wall', P['fins']: 'wall'}
    types = {k: v for k, v in types.items() if k in quads and len(quads[k])}
    walls = [k for k, v in types.items() if v == 'wall']
    info = dict(sector=MP.SECTOR, symmetryFraction=frac, fins=bool(MP.FINS_ON),
                nCells=int(n_cells), hScale=float(MP.H_SCALE),
                yPlusTarget=float(MP.YPLUS_TARGET), y1=float(d['y1']),
                Uref=float(MP.U), nuRef=float(MP.NU),
                Aref=float(ref['Aref']), lRef=float(ref['lRef']),
                rBody=float(G.R_BODY), lBody=float(G.L_TOTAL))
    # Fin bounding geometry, so case-*/system/topoSetDict can place a local
    # refinement region without repeating any of these numbers.  See
    # case-*/Allrefine and docs/WORKFLOW.md section 3.
    if MP.FINS_ON:
        info.update(finRootR=float(G.FIN_ROOT_R), finTipR=float(G.FIN_TIP_R),
                    finTipSmear=float(MP.FIN_TIP_SMEAR),
                    finX0=float(G.FIN_ROOT_LE), finX1=float(G.FIN_TIP_TE))
    info.update(wallPatches=walls, patchTypes=types)
    return info


def main(argv=None):
    a = parse_args(argv)
    prov = apply_cli(a)
    verbose = not a.quiet

    if a.plan:
        MP.report()
        return 0

    import meshFinish
    import sectorAssembly
    import meshIO
    import finPatch

    t0 = time.time()
    d = MP.derived()                              # validates the parameter set
    if verbose:
        print(f'building  sector={MP.SECTOR}  H_SCALE={MP.H_SCALE:g}  fins={MP.FINS_ON}')
        for p in prov:
            print('  ' + p)
        est = MP.predicted_cells(d)
        print(f"  predicted cells: {est['quadrant']:,d} per quadrant, {est['total']:,d} total")

    q = meshFinish.generate_quadrant(verbose)
    P = MP.PATCHES
    n_copies = MP.n_copies()

    # the quadrant is always audited: cheap, and every copy has the same quality
    symm_q, fin_q = finPatch.split_symm(q['quads'][P['symmetry']], q['fin_nodes'])
    quads_q = dict(q['quads'])
    quads_q[P['symmetry']] = symm_q
    quads_q[P['fins']] = fin_q
    qa = meshFinish.audit(q['nodes'], q['hexes'], quads_q, verbose, label='quadrant audit')
    if MP.FINS_ON and verbose:
        a_fin = finPatch.wetted_area(q['nodes'], fin_q) / 2.0      # two half fins
        a_ref = finPatch.wetted_area_analytic()                    # the exact surface
        print(f'  fin: {len(fin_q):,d} faces, wetted area per half-fin '
              f'{a_fin*1e4:.2f} cm2 vs {a_ref*1e4:.2f} exact '
              f'({100.0*(a_fin-a_ref)/a_ref:+.1f} %)')
        if MP.fin_edge_fit():
            print(f'       edges fitted: the outline IS the planform.  What is '
                  f'left is the reference, which integrates from FIN_ROOT_R,')
            print(f'       0.5 mm inside the wall the mesh starts at, and on '
                  f'y = r rather than the cylinder the fin is wrapped on')
        else:
            print(f'       the excess is the one-cell rim where the patch outline '
                  f'is quantised -- set FIN_EDGE_FIT to put the edges on the mesh')

    if n_copies == 1:
        nodes, hexes, quads = q['nodes'], q['hexes'], quads_q
        qf = qa
    else:
        nodes, hexes, quads = sectorAssembly.assemble(
            q['nodes'], q['hexes'], q['quads'], q['plane_a'], q['plane_b'],
            q['fin_nodes'], n_copies, verbose=verbose)
        if a.no_audit:
            qf = dict(qa, n_cells=len(hexes), audited=False)
            if verbose:
                print('  (assembled mesh not audited: --no-audit)')
        else:
            qf = meshFinish.audit(nodes, hexes, quads, verbose,
                                  label=f'{MP.SECTOR} audit')

    if not qf['patches_consistent']:
        print('!! boundary faces and patch faces disagree -- refusing to write', file=sys.stderr)
        return 2
    if qf['n_negative'] or qf['triple_face']:
        print('!! degenerate cells or faces -- refusing to write', file=sys.stderr)
        return 2

    if a.no_write:
        if verbose:
            print(f'  nothing written (--no-write).  total {time.time()-t0:.1f}s')
        return 0

    name = a.name or auto_name()
    out = a.out or os.path.join(HERE, 'output', name + '.msh')
    base = out[:-4] if out.endswith('.msh') else out
    os.makedirs(os.path.dirname(os.path.abspath(out)) or '.', exist_ok=True)

    order = [P['inlet'], P['outlet'], P['farfield'], P['symmetry'],
             P['nose'], P['body'], P['tail'], P['fins']]
    order = [k for k in order if k in quads] + [k for k in quads if k not in order]
    tw = time.time()
    meshIO.write_msh2(out, nodes, hexes, quads, order)
    info = mesh_info(len(hexes), d, quads)
    meshIO.write_meshinfo(base + '.meshInfo', info)
    meshIO.write_params(base + '.params.py', MP.snapshot(),
                        header=f'{name}: {len(hexes):,d} cells\n' + '\n'.join(prov))
    if verbose:
        print(f'\n  wrote {out}  ({os.path.getsize(out)/1e6:.0f} MB, v2.2 ASCII, '
              f'{time.time()-tw:.0f}s)')
        print(f'        {base}.meshInfo   (Aref {info["Aref"]:.6g} m2, '
              f'patches {", ".join(k for k in order if len(quads[k]))})')
        print(f'        {base}.params.py  (rebuild with --preset)')
        print(f'  total {time.time()-t0:.1f}s')
    return 0


if __name__ == '__main__':
    sys.exit(main())
