# case-central/ — case template, M > 1

`rhoCentralFoam` with the Kurganov central-upwind flux, `kOmegaSST` with wall
functions, and air as a perfect gas with Sutherland viscosity. The sweep uses
it for M 1.2 and 1.8, started from a `case-pimple/` run at the same flight
condition. M < 1 uses `case-pimple/`, which takes its shared files from here.

> **CAUTION** The template is not validated, and the central run at M 1.8 does
> not converge yet (user guide, Troubleshooting). Do not use the results before
> you do the procedures of the chapter Validation.

Do not run a case in this directory. Make a case with `../newCase.sh`:

```bash
../newCase.sh ~/runs/m18 --mesh ~/meshes/h1 --Minf 1.8 --init ~/runs/m18-pimple
```

| Item | Function | User guide |
|---|---|---|
| `system/flowConditions` | The only file that sets the flight condition. | Solvers and models > Flight condition |
| `config.json` | Courant ramp and stop criterion. `newCase.sh` renders the `*.j2` files with it. | Solvers and models > config.json |
| `0.orig/p` | The only boundary condition that changes with the Mach number, at `Minf = 1`. | Solvers and models > Boundary conditions |
| `system/FOs/FOcommon` | Function objects of the two templates, including the samples for Nektar++. | Results and measurement |
| `constant/meshInfo` | `newCase.sh` replaces it with the `meshInfo` of the mesh. | Mesh > meshInfo |

The user guide is `../docs/guide/aconcagua-user-guide.pdf`.
