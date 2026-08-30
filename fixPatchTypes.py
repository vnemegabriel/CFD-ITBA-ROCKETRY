#!/usr/bin/env python3
"""
fixPatchTypes.py -- run this IMMEDIATELY after gmshToFoam. Not optional.

gmshToFoam creates every patch as plain `type patch`, because a .msh physical
group carries a NAME but no notion of what kind of boundary it is.  That is a
silent, load-bearing problem:

  * `symm` must be `type symmetry`, or `0/U`'s `type symmetry;` entry raises
    a type mismatch and the solver refuses to start.  (That one at least fails
    loudly.)
  * `cone`, `walls` and `tail` must be `type wall`, and THIS one does not fail
    loudly.  nutUSpaldingWallFunction on a non-wall patch, wall distance, yPlus
    and forceCoeffs all quietly do the wrong thing.  You would get a converged
    run and a wrong Cd.

So the conversion is gmshToFoam followed by this, always, as one step.
"""

import re
import sys

TYPES = {
    'cone':   'wall',            # nose
    'walls':  'wall',            # cylindrical body
    'tail':   'wall',            # boattail + flat base
    'fins':   'wall',            # the two half fins in the symmetry planes
    'symm':   'symmetry',        # the two quarter-symmetry planes
    'inlet':  'patch',
    'outlet': 'patch',
    'box':    'patch',           # cylindrical farfield (name kept for parity)
}

# The boundary file opens with a FoamFile header that looks exactly like a
# patch block to a regex.  It is not one.
SKIP = {'FoamFile'}

BLOCK = re.compile(r'^(\s*)(\w+)\s*\n(\s*)\{(.*?)^\3\}', re.S | re.M)


def main(path='constant/polyMesh/boundary'):
    src = open(path).read()
    seen, unknown, changed = [], [], 0

    def fix(m):
        nonlocal changed
        ind, name, bind, body = m.groups()
        if name in SKIP:
            return m.group(0)
        if name not in TYPES:
            unknown.append(name)
            return m.group(0)
        seen.append(name)
        want = TYPES[name]
        # gmshToFoam also writes a legacy `physicalType`; drop it rather than
        # leave a second, contradictory type declaration in the file
        body = re.sub(r'\n\s*physicalType[^\n]*\n', '\n', body)
        new, n = re.subn(r'(\btype\s+)\w+(\s*;)', rf'\g<1>{want}\g<2>', body, count=1)
        if n == 0:
            new = f'\n{bind}    type            {want};' + body
        if want == 'wall':
            new = re.sub(r'\n\s*inGroups[^\n]*\n', '\n', new)
            new = new.replace(f'type            {want};',
                              f'type            {want};\n{bind}    inGroups        1(wall);', 1)
        if new != body:
            changed += 1
        return f'{ind}{name}\n{bind}{{{new}{bind}}}'

    out = BLOCK.sub(fix, src)
    open(path, 'w').write(out)

    missing = sorted(set(TYPES) - set(seen))
    print(f'  patch types set on {changed} of {len(seen)} patches')
    for n in seen:
        print(f'    {n:<10} -> {TYPES[n]}')
    if unknown:
        print(f'  !! UNRECOGNISED PATCHES: {unknown}')
        if 'defaultFaces' in unknown:
            print('     defaultFaces means some boundary faces had no physical group')
            print('     in the .msh. Do not proceed -- fix the export.')
        else:
            print('     Not defaultFaces, so the export is fine -- this script is')
            print('     just missing an entry. Add it to TYPES above.')
        return 1
    if missing:
        print(f'  !! EXPECTED BUT ABSENT: {missing}')
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:]))
