# case-pimple/ — case template, M < 1 and the start of M > 1

`rhoPimpleFoam` with `transonic yes` and local time stepping (`localEuler`).
Same turbulence model, gas, boundary conditions, function objects and flight
condition as `case-central/`. For M > 1 it makes the initial field of the
central case: from the uniform free stream the central case stops at the base
(user guide, Solvers and models > Start of the supersonic cases).

This directory holds only the files that differ from `case-central/`.
`../newCase.sh --solver pimple` copies `case-central/` without the files in
`central.exclude`, and then copies this directory on top.

```bash
../newCase.sh ~/runs/m03-h1 --mesh ~/meshes/h1 --Minf 0.3 --solver pimple
../newCase.sh ~/runs/m03-h2 --mesh ~/meshes/h2 --Minf 0.3 --solver pimple --init ~/runs/m03-h1
```

| File | Function | User guide |
|---|---|---|
| `config.json` | `maxCo`, `minIter`, `endIter`, stop criterion in iterations | Solvers and models > config.json |
| `system/fvSolution.j2` | `PIMPLE`, including the local-time-step controls | rhoPimpleFoam template > PIMPLE |
| `system/fvSchemes` | `localEuler` and the bounded second-order schemes | rhoPimpleFoam template > Schemes |
| `system/fvOptions` | `limitTemperature`, 100–1000 K | rhoPimpleFoam template |
| `system/controlDict.j2` | Time in iterations; `FOs/FOcommon` and `FOs/FOrunTimeControl.j2` from `case-central/` | Solvers and models > Local time step |
| `central.exclude` | Files of `case-central/` that this template replaces or does not use | |

To change a file that both templates use (`0.orig/`, `constant/`,
`flowConditions`, `FOs/FOcommon`), change it in `case-central/`.

The user guide is `../docs/guide/aconcagua-user-guide.pdf`.
