#!/usr/bin/env python3
"""
fixPatchTypes.py -- run IMMEDIATELY after gmshToFoam.  Allmesh does; it is
not optional.

gmshToFoam creates every patch as plain `type patch`, because a .msh physical
group carries a NAME but no notion of what kind of boundary it is.  That is a
silent, load-bearing problem:

  * `symm` must be `type symmetry`, or `0/U`'s `type symmetry;` entry raises a
    type mismatch and the solver refuses to start.  (That one fails loudly.)
  * `cone`, `walls`, `tail`, `fins` must be `type wall`, and THIS does not fail
    loudly.  nutUSpaldingWallFunction on a non-wall patch, wall distance, yPlus
    and forceCoeffs all quietly do the wrong thing.  You would get a converged
    run and a wrong Cd.

The types come from constant/meshInfo (the `patchTypes` dictionary that
mesh/build.py writes for every mesh), so a mesh without a `symm` patch -- a
full 360 deg one -- is handled without editing anything here.  The hard-coded
table below is only the fallback for a mesh that arrived without its sidecar.

Exit status is non-zero on anything unexpected: a patch the mesh has that the
table does not, or the other way round.  Allmesh stops on that.
"""

import re
import sys

FALLBACK = {
    'cone': 'wall', 'walls': 'wall', 'tail': 'wall', 'fins': 'wall',
    'symm': 'symmetry', 'inlet': 'patch', 'outlet': 'patch', 'box': 'patch',
}

SKIP = {'FoamFile'}     # the header looks like a patch block to a regex; it is not

BLOCK = re.compile(r'^(\s*)(\w+)\s*\n(\s*)\{(.*?)^\3\}', re.S | re.M)


def read_types(info='constant/meshInfo'):
    """The patchTypes sub-dictionary of meshInfo, or the fallback table."""
    try:
        src = open(info).read()
    except OSError:
        print(f'  (no {info}: using the built-in patch table)')
        return dict(FALLBACK), False
    src = re.sub(r'/\*.*?\*/', '', src, flags=re.S)
    src = re.sub(r'//[^\n]*', '', src)
    m = re.search(r'patchTypes\s*\{(.*?)\}', src, re.S)
    if not m:
        print(f'  ({info} has no patchTypes entry: using the built-in patch table)')
        return dict(FALLBACK), False
    types = dict(re.findall(r'(\w+)\s+(\w+)\s*;', m.group(1)))
    if not types:
        raise SystemExit(f'  !! patchTypes in {info} is empty')
    return types, True


def main(path='constant/polyMesh/boundary', info='constant/meshInfo'):
    types, from_info = read_types(info)
    src = open(path).read()
    seen, unknown, changed = [], [], 0

    def fix(m):
        nonlocal changed
        ind, name, bind, body = m.groups()
        if name in SKIP:
            return m.group(0)
        if name not in types:
            unknown.append(name)
            return m.group(0)
        seen.append(name)
        want = types[name]
        # gmshToFoam also writes a legacy `physicalType`; drop it rather than
        # leave a second, contradictory type declaration in the file
        body = re.sub(r'\n\s*physicalType[^\n]*\n', '\n', body)
        new, n = re.subn(r'(\btype\s+)\w+(\s*;)', rf'\g<1>{want}\g<2>', body, count=1)
        if n == 0:
            new = f'\n{bind}    type            {want};' + body
        new = re.sub(r'\n\s*inGroups[^\n]*\n', '\n', new)
        if want == 'wall':
            new = new.replace(f'type            {want};',
                              f'type            {want};\n{bind}    inGroups        1(wall);', 1)
        if new != body:
            changed += 1
        return f'{ind}{name}\n{bind}{{{new}{bind}}}'

    out = BLOCK.sub(fix, src)
    open(path, 'w').write(out)

    missing = sorted(set(types) - set(seen))
    print(f'  patch types set on {changed} of {len(seen)} patches'
          f' ({"from " + info if from_info else "built-in table"})')
    for n in seen:
        print(f'    {n:<10} -> {types[n]}')
    if unknown:
        print(f'  !! UNRECOGNISED PATCHES: {unknown}')
        if 'defaultFaces' in unknown:
            print('     defaultFaces means some boundary faces had no physical group')
            print('     in the .msh.  Do not proceed -- rebuild the mesh.')
        else:
            print('     The mesh has a patch its meshInfo does not list; rebuild the')
            print('     mesh so the sidecar matches, or add it to FALLBACK here.')
        return 1
    if missing:
        print(f'  !! EXPECTED BUT ABSENT: {missing}')
        print('     meshInfo lists a patch the mesh does not have -- wrong sidecar?')
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))
