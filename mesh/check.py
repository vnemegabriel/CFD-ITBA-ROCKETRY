#!/usr/bin/env python3
"""
check.py -- the regression battery.  Run it before saying a change is done.

    python check.py             fast   (~1 min)  smoke + coarse quadrant
    python check.py --full      full   (~6 min)  + sector assembly, --no-fins,
                                                 every preset, the .params.py
                                                 round-trip, and `fine`
    python check.py --update    run and REWRITE the baseline
    python check.py --run NAME  internal: build one config, print its metrics

WHAT IT CHECKS, AND WHY IN TWO LAYERS
=====================================
An absolute FLOOR that never moves, because these are the things that make a
mesh unusable rather than merely worse:

    cells with volume <= 0        a face shared by 3+ cells
    boundary faces == patch faces (nothing in defaultFaces, every stitched
                                   face really became interior)

and a BASELINE, check.baseline.json, holding the numbers this repo currently
produces.  A fixed quality threshold would not have caught the two regressions
that actually happened while this mesh was being built: the fin patch coming
out 10-21 % larger than the planform, and the mean non-orthogonality moving
from 5.78 to 7.59 deg.  Neither breaks a threshold; both are visible instantly
against the previous run.

So: the floor says "this mesh is broken", the baseline says "this mesh is
different, come and explain why".  When the difference is intended -- it
usually is -- `--update` records the new numbers and the diff shows up in the
commit next to the change that caused it.

CONFIGURATIONS
==============
Chosen to exercise the paths that broke at least once each:

    smoke           the pipeline, in seconds
    coarse          the working mesh
    defaults        meshParams as it stands, --plan only: catches a default
                    left in a state that cannot even be validated
    coarse+fintip   the tip zone: the node line at FIN_TIP_R
    half            sectorAssembly: nodes merged, faces made interior, the fin
                    faces NOT merged
    nofins          the FINS_ON = False path, which no preset covers
    norelax         AZ_RELAX = False.  Freezes the mesh as it was before the
                    azimuthal relaxation existed, so the relaxation can always
                    be told apart from whatever else moved
    presets         every preset still validates
    roundtrip       a .params.py sidecar rebuilds the mesh it describes
    fine            the reference resolution
"""

import argparse
import json
import os
import runpy
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
BASELINE = os.path.join(HERE, 'check.baseline.json')

# ---------------------------------------------------------------- configs --
# (name, [presets], {overrides}, sector)
BUILDS = {
    'smoke':         ([], dict(H_SCALE=6.0, FIN_H_R=None), 'quarter'),
    'coarse':        (['coarse'], {}, 'quarter'),
    'coarse+fintip': (['coarse', 'fintip'], {}, 'quarter'),
    'half':          (['coarse', 'fintip'], {}, 'half'),
    'nofins':        (['coarse'], dict(FINS_ON=False), 'quarter'),
    'norelax':       (['coarse'], dict(AZ_RELAX=False), 'quarter'),
    'fine':          (['fine', 'fintip'], {}, 'quarter'),
}
FAST = ['smoke', 'coarse']
FULL = ['smoke', 'coarse', 'coarse+fintip', 'half', 'nofins', 'norelax', 'fine']

PRESETS = ['coarse', 'medium', 'fine', 'fintip', 'wake_unsteady', 'smoke']

# metric -> (rule, tolerance).  'up' fails only when the number GROWS, 'both'
# when it moves either way.  Chosen so ordinary re-solving noise passes and a
# real change does not.
TOL = {
    'n_cells':          ('both', 0.005),
    'vol_total':        ('both', 1e-6),   # the domain itself must not move
    'nonortho_max':     ('up',   0.01),
    'nonortho_mean':    ('up',   0.10),
    'nonortho_gt70':    ('up',   0.05),
    'nonortho_gt40':    ('up',   0.25),   # the fin-edge shear moves this a lot
    'skew_max':         ('up',   0.02),
    'ar_max':           ('up',   0.05),
    'fin_faces':        ('both', 0.02),
}
# percentage POINTS, not relative: the fin area error is already a percentage
TOL_ABS = {'fin_area_err_pct': 0.30}


# ------------------------------------------------------------- one config --
def run_one(name):
    """Build one configuration and return its metrics.  Runs in a child."""
    sys.path.insert(0, HERE)
    import meshParams as MP
    presets, over, sector = BUILDS[name]
    for p in presets:
        ns = runpy.run_path(os.path.join(HERE, 'presets', p + '.py'))
        MP.apply_overrides({k: v for k, v in ns.items()
                            if k.isupper() and not k.startswith('_')})
    if over:
        MP.apply_overrides(over)
    MP.apply_overrides(dict(SECTOR=sector))

    import meshFinish, finPatch, sectorAssembly
    t0 = time.time()
    q = meshFinish.generate_quadrant(verbose=False)
    P = MP.PATCHES
    symm, fin = finPatch.split_symm(q['quads'][P['symmetry']], q['fin_nodes'])
    quads = dict(q['quads'])
    quads[P['symmetry']], quads[P['fins']] = symm, fin

    if MP.n_copies() == 1:
        nodes, hexes = q['nodes'], q['hexes']
    else:
        nodes, hexes, quads = sectorAssembly.assemble(
            q['nodes'], q['hexes'], q['quads'], q['plane_a'], q['plane_b'],
            q['fin_nodes'], MP.n_copies(), verbose=False)
    a = meshFinish.audit(nodes, hexes, quads, verbose=False)

    m = {k: a[k] for k in ('n_cells', 'n_internal', 'n_boundary', 'n_negative',
                           'vol_total', 'nonortho_max', 'nonortho_mean',
                           'nonortho_gt70', 'nonortho_gt40', 'skew_max',
                           'ar_max')}
    m['triple_face'] = bool(a['triple_face'])
    m['patches_consistent'] = bool(a['patches_consistent'])
    m['fin_faces'] = len(quads[P['fins']])
    if MP.FINS_ON:
        got = finPatch.wetted_area(nodes, quads[P['fins']]) / (2.0 * MP.n_copies())
        ref = finPatch.wetted_area_analytic()
        m['fin_area_err_pct'] = 100.0 * (got - ref) / ref
    m['seconds'] = round(time.time() - t0, 1)
    return m


# ------------------------------------------------------------- comparison --
def floor_failures(name, m):
    """The things that make a mesh unusable.  Independent of any baseline."""
    bad = []
    if m['n_negative']:
        bad.append(f'{m["n_negative"]:,d} cells with volume <= 0')
    if m['triple_face']:
        bad.append('a face is shared by 3 or more cells')
    if not m['patches_consistent']:
        bad.append('boundary faces and patch faces disagree (defaultFaces)')
    if m['n_boundary'] == 0 or m['n_cells'] == 0:
        bad.append('empty mesh')
    return bad


def drifts(m, base):
    """Metrics that moved further than TOL.  Returns a list of strings."""
    out = []
    for k, (rule, tol) in TOL.items():
        if k not in m or k not in base:
            continue
        a, b = float(base[k]), float(m[k])
        if abs(a) < 1e-30:
            moved = abs(b) > 1e-30
            rel = float('inf') if moved else 0.0
        else:
            rel = (b - a) / abs(a)
        if (rule == 'up' and rel > tol) or (rule == 'both' and abs(rel) > tol):
            out.append(f'{k}: {a:,.4g} -> {b:,.4g}  ({rel*100:+.1f} %)')
    for k, tol in TOL_ABS.items():
        if k in m and k in base and abs(float(m[k]) - float(base[k])) > tol:
            out.append(f'{k}: {float(base[k]):+.2f} -> {float(m[k]):+.2f} '
                       f'(> {tol} points)')
    return out


# ------------------------------------------------------------ subprocesses --
def build_in_child(name):
    r = subprocess.run([sys.executable, os.path.abspath(__file__), '--run', name],
                       capture_output=True, text=True, cwd=HERE)
    if r.returncode != 0:
        return None, (r.stderr.strip().splitlines() or ['(no output)'])[-1]
    try:
        return json.loads(r.stdout.strip().splitlines()[-1]), None
    except (ValueError, IndexError):
        return None, 'could not parse the child output'


def plan_only(args, label):
    r = subprocess.run([sys.executable, os.path.join(HERE, 'build.py'), *args, '--plan'],
                       capture_output=True, text=True, cwd=HERE)
    if r.returncode == 0:
        return True, ''
    tail = (r.stderr.strip().splitlines() or [''])[-1]
    return False, tail[:150]


def roundtrip():
    """A .params.py sidecar must rebuild the mesh it describes."""
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, 'rt.msh')
        r = subprocess.run([sys.executable, os.path.join(HERE, 'build.py'),
                            '--preset', 'coarse', '--out', out, '--quiet'],
                           capture_output=True, text=True, cwd=HERE)
        if r.returncode != 0:
            return False, 'the first build failed'
        side = out[:-4] + '.params.py'
        if not os.path.exists(side):
            return False, 'no .params.py was written'
        r = subprocess.run([sys.executable, os.path.join(HERE, 'build.py'),
                            '--preset', side, '--plan'],
                           capture_output=True, text=True, cwd=HERE)
        if r.returncode != 0:
            return False, 'the sidecar does not rebuild: ' + \
                (r.stderr.strip().splitlines() or [''])[-1][:120]
    return True, ''


# -------------------------------------------------------------------- main --
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--full', action='store_true', help='the whole battery')
    ap.add_argument('--update', action='store_true',
                    help='rewrite check.baseline.json from this run')
    ap.add_argument('--run', metavar='NAME', help=argparse.SUPPRESS)
    a = ap.parse_args(argv)

    if a.run:                                   # child process
        print(json.dumps(run_one(a.run)))
        return 0

    base = {}
    if os.path.exists(BASELINE):
        base = json.load(open(BASELINE))
    elif not a.update:
        print(f'no baseline at {BASELINE}: run  python check.py --full --update')

    names = FULL if a.full else FAST
    t0 = time.time()
    results, failed, drifted = {}, [], []

    # flush as we go: --full takes minutes and a piped run would otherwise
    # show nothing at all until it finished
    print(f'{"config":<16}{"cells":>11}{"nOrt max":>10}{"mean":>7}{"skew":>7}'
          f'{"fin area":>10}{"s":>6}  status', flush=True)
    print('-' * 78, flush=True)
    for n in names:
        m, err = build_in_child(n)
        if m is None:
            print(f'{n:<16}{"--":>11}{"":>40}  BUILD FAILED: {err[:40]}',
                  flush=True)
            failed.append((n, [err]))
            continue
        results[n] = m
        bad = floor_failures(n, m)
        d = drifts(m, base.get(n, {})) if n in base else []
        area = f"{m['fin_area_err_pct']:+.2f} %" if 'fin_area_err_pct' in m else '--'
        status = 'FLOOR' if bad else ('drift' if d else 'ok')
        print(f'{n:<16}{m["n_cells"]:>11,d}{m["nonortho_max"]:>10.2f}'
              f'{m["nonortho_mean"]:>7.2f}{m["skew_max"]:>7.3f}{area:>10}'
              f'{m["seconds"]:>6.0f}  {status}', flush=True)
        if bad:
            failed.append((n, bad))
        if d:
            drifted.append((n, d))

    if a.full:
        print()
        for p in PRESETS:
            ok, why = plan_only(['--preset', p], p)
            print(f'  preset {p:<16}{"ok" if ok else "FAILED: " + why}',
                  flush=True)
            if not ok:
                failed.append((f'preset {p}', [why]))
        ok, why = plan_only([], 'defaults')
        print(f'  meshParams defaults    {"ok" if ok else "FAILED: " + why}')
        if not ok:
            failed.append(('defaults', [why]))
        ok, why = roundtrip()
        print(f'  .params.py round-trip  {"ok" if ok else "FAILED: " + why}')
        if not ok:
            failed.append(('roundtrip', [why]))

    if a.update and failed:
        # Never bake a broken mesh into the reference: the whole point of the
        # baseline is that it describes a mesh someone was willing to ship.
        print('\nNOT writing the baseline: the run has failures (see below).')
    elif a.update:
        merged = dict(base)
        merged.update(results)
        with open(BASELINE, 'w', newline='\n') as f:
            json.dump(merged, f, indent=2, sort_keys=True)
            f.write('\n')
        print(f'\nbaseline written: {BASELINE}')

    print(f'\ntotal {time.time()-t0:.0f}s')
    if drifted and not a.update:
        print('\nMOVED against the baseline.  Intended?  Then rerun with '
              '--update and commit the baseline with the change:')
        for n, d in drifted:
            for line in d:
                print(f'  {n:<16}{line}')
    if failed:
        print('\nFAILED:')
        for n, why in failed:
            for line in why:
                print(f'  {n:<16}{line}')
        return 2
    if drifted and not a.update:
        return 1
    print('\nall good.' if not a.update else '')
    return 0


if __name__ == '__main__':
    sys.exit(main())
