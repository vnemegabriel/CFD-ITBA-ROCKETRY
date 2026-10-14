# case-central/ — case template, M > 1

`rhoCentralFoam` with the Kurganov central-upwind flux, `kOmegaSST` with wall
functions, and air as a perfect gas with Sutherland viscosity. The sweep uses
it for M 1.2 and 1.8, started from a `case-pimple/` run at the same Mach
(docs/SOLVERS.md §17). M < 1 uses `case-pimple/`, which takes its shared
files from here.

> **CAUTION** The template is not validated. Do not use the results before
> you do the procedures in [../docs/VALIDATION.md](../docs/VALIDATION.md).

Do not run a case in this directory. Make a case with `../newCase.sh`:

```bash
../newCase.sh ~/runs/m18 --mesh ~/meshes/base --Minf 1.8
```

| Item | Function | Reference |
|---|---|---|
| `system/flowConditions` | The only file that sets the flight condition. | [SOLVERS.md §3](../docs/SOLVERS.md) |
| `config.json` | Courant ramp and stop criterion. `newCase.sh` renders the `*.j2` files with it. | [SOLVERS.md §5](../docs/SOLVERS.md) |
| `0.orig/p` | The only boundary condition that changes with the Mach number. It changes at `Minf = 1`. | [SOLVERS.md §9](../docs/SOLVERS.md) |
| `constant/meshInfo` | `newCase.sh` replaces it with the `meshInfo` of the mesh. | [MESH.md §7](../docs/MESH.md) |

The mesh does not set y+, and the velocity changes by a factor of 6 in the
sweep. Measure y+ on each case.
