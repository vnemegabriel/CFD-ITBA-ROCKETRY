# Workflow

This document gives the procedures that change the STL geometry into force
coefficients. Each procedure gives the command, the steps that the script
does, and the result.

- [MESH.md](MESH.md) gives the mesh dictionaries.
- [SOLVERS.md](SOLVERS.md) gives the case dictionaries and the models.
- [VALIDATION.md](VALIDATION.md) gives the validation procedures.

## 1. Directory structure

| Path | Contents |
|---|---|
| `mesh/Allmesh` | Script that makes one mesh. |
| `mesh/stl/Aconcagua.stl` | Geometry of the rocket. |
| `mesh/system/` | `meshDict`, `createPatchDict`, and the `controlDict`, `fvSchemes`, `fvSolution` that the mesh applications need. |
| `mesh/meshInfo` | Reference values and fin geometry of the mesh. |
| `case-central/` | Case template for `rhoCentralFoam`. |
| `common/Allrun`, `common/Allrefine` | Scripts that `newCase.sh` copies into each case. |
| `newCase.sh` | Script that makes one case from the template and one mesh. |
| `run.sh`, `sweep.txt` | Script and input file for a sweep of cases. |

## 2. Requirements

| Item | Requirement |
|---|---|
| OpenFOAM | ESI v2412, with cfMesh (`cartesianMesh`). |
| `jinja2` | `jinja2-cli`. Install it with `pipx install jinja2-cli`. |
| Shell tools | `bash`, `rsync`, `awk`, GNU `sed`. |

Do this step before you use a script:

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
```

All scripts stop with the message `source the OpenFOAM bashrc first` if you
do not do this step.

> **CAUTION** Put the mesh directories and the case directories outside
> OneDrive and outside `/mnt/c`. These file systems make the OpenFOAM input
> and output several times slower.

## 3. Make a mesh

### 3.1 Command

```bash
mesh/Allmesh <mesh-dir> [stl-file] [--maxCellSize <m>] [--finLevel <n>]
```

| Argument | Default | Function |
|---|---|---|
| `<mesh-dir>` | none, mandatory | Directory for the new mesh. It must not exist. |
| `[stl-file]` | `mesh/stl/Aconcagua.stl` | Geometry. MESH.md §3 gives the requirements. |
| `--maxCellSize <m>` | value in `meshDict` (1.6 m) | Sets `maxCellSize` in the copy of `meshDict`. |
| `--finLevel <n>` | value in `meshDict` (9) | Sets the refinement level of the `fins` patch in the copy of `meshDict`. |

### 3.2 Procedure

1. Run `mesh/Allmesh <mesh-dir>`.
2. Read the summary at the end of the output. It gives the number of cells,
   the cell types, the maximum non-orthogonality, the maximum skewness and
   the failed checks.
3. If the output shows `negative cell volumes`, do not use the mesh. Change
   `meshDict` and make a new mesh in a new directory.
4. Open `<mesh-dir>/constant/polyMesh/boundary`. Make sure that the patches
   `cone`, `walls`, `tail` and `fins` have the type `wall`. Refer to MESH.md
   §6.

### 3.3 Result

| Item | Contents |
|---|---|
| `<mesh-dir>/constant/polyMesh/` | The mesh. |
| `<mesh-dir>/system/` | A copy of `mesh/system/`, with the options applied. Each mesh keeps the `meshDict` that made it. |
| `<mesh-dir>/meshInfo` | A copy of `mesh/meshInfo`. `Allmesh` adds `nCells`, `maxCellSize`, `finLevel` and `meshDir`. |
| `<mesh-dir>/log.*` | One log for each application. |

`Allmesh` does not write into a directory that exists. To make the mesh
again, remove the directory or use a new name.

## 4. Make a case

### 4.1 Command

```bash
./newCase.sh <run-dir> --mesh <mesh-dir> [options]
```

| Option | Function |
|---|---|
| `--mesh <mesh-dir>` | Mandatory. The path of a mesh that `Allmesh` made. |
| `--Minf <M>` | Free-stream Mach number. |
| `--pInf <Pa>` | Free-stream static pressure. |
| `--Tinf <K>` | Free-stream static temperature. |
| `--alpha <deg>` | Angle of attack. |
| `--beta <deg>` | Angle of sideslip. |
| `--Ti <->` | Free-stream turbulence intensity. |
| `--nuRatio <->` | Ratio `nut/nu` in the free stream. |
| `--np <N>` | Number of subdomains in `system/decomposeParDict`. |
| `--refine tip` or `--refine fins` | Runs `./Allrefine` on the new case. Refer to §7. |

The flow options change entries in `system/flowConditions`. The options that
you do not give keep the template values.

### 4.2 Steps that newCase.sh does

1. It stops if `--mesh` is missing, if the mesh directory has no
   `constant/polyMesh`, or if `<run-dir>` exists.
2. It copies `case-central/` into `<run-dir>`. It does not copy the mesh,
   the logs, `0/`, `processor*`, `postProcessing/` or the time directories.
3. It copies `<mesh-dir>/constant/polyMesh/` into `<run-dir>/constant/`.
4. It copies `<mesh-dir>/meshInfo` to `<run-dir>/constant/meshInfo`. This
   file replaces the `meshInfo` of the template.
5. It copies `common/Allrun` and `common/Allrefine` into `<run-dir>`.
6. It renders each `*.j2` file with the values in `config.json`, and then
   removes the `*.j2` file. Refer to SOLVERS.md §5.
7. It writes each flow option into `system/flowConditions`. It stops if the
   entry does not exist in the file.
8. It writes `--np` into `system/decomposeParDict`.
9. It runs `./Allrefine` if you give `--refine`.

> **NOTE** In `newCase.sh`, `--mesh` is a path. Only `run.sh` accepts a mesh
> name without a path.

## 5. Run a case

### 5.1 Command

```bash
cd <run-dir>
./Allrun        # parallel, numberOfSubdomains from system/decomposeParDict
./Allrun 16     # parallel on 16 processes; writes 16 into decomposeParDict
./Allrun 1      # serial
```

### 5.2 Steps that Allrun does

1. It stops if `constant/polyMesh` does not exist.
2. It compares `alpha` and `beta` in `system/flowConditions` with `sector`
   in `constant/meshInfo`. It stops if the mesh cannot show the flow angle.
3. It copies `0.orig/` to `0/` (`restore0Dir`).
4. In parallel: it runs `decomposePar -force`, then `rhoCentralFoam`, then
   `reconstructPar -latestTime`. In serial: it runs `rhoCentralFoam`.
5. It shows the last line of the force coefficients.

| `sector` in `meshInfo` | Permitted angles |
|---|---|
| `quarter` | `alpha = 0` and `beta = 0` |
| `half` | `beta = 0` |
| `full` | All angles |

The mesh that `Allmesh` makes is `full`.

### 5.3 End of a run

The run stops at `endTime` or when the `runTimeControl` function object finds
that the average of `Cd` is constant. SOLVERS.md §5 gives the criterion and
its limits.

### 5.4 Continue a run

`controlDict` has `startFrom latestTime`. Run `./Allrun` again in the same
directory to continue from the last time directory.

> **CAUTION** `Allrun` runs `decomposePar -force`, which removes the
> `processor*` directories. If the solver stopped before `reconstructPar`,
> run `reconstructPar -latestTime` before you run `./Allrun` again. If you do
> not, the data after the last reconstructed time is lost.

On the first continuation, read the start time in `log.rhoCentralFoam`. Make
sure that it is the last time and not 0.

### 5.5 Start a run again from zero

1. Remove the time directories and the processor directories:
   `rm -rf [1-9]* processor*`.
2. Run `./Allrun`.

`./Allclean` also removes `constant/polyMesh`. After `./Allclean`, make a new
case with `newCase.sh`.

## 6. Run a sweep

### 6.1 Input file

`sweep.txt` has one case on each line: a case name, then the options of
`newCase.sh`.

```
m03     --mesh base --Minf 0.3
m08     --mesh base --Minf 0.8
m03-h2  --mesh h2 --maxCellSize 1.13 --Minf 0.3
```

- Lines that start with `#`, and empty lines, have no effect.
- `--mesh` is mandatory on each line.
- A mesh name without `/` is a directory in `<meshes-dir>`.
- `--maxCellSize` and `--finLevel` apply only when `run.sh` makes the mesh.

### 6.2 Command

```bash
./run.sh [sweep-file] [runs-dir] [meshes-dir]
```

| Argument | Default |
|---|---|
| `sweep-file` | `flightModel/sweep.txt` |
| `runs-dir` | `flightModel/runs` |
| `meshes-dir` | `~/meshes` |

> **CAUTION** The default `runs-dir` is in the repository. Give a
> `runs-dir` outside the repository, outside OneDrive and outside `/mnt/c`.

### 6.3 Steps that run.sh does for each line

1. If `<runs-dir>/<name>` exists, it goes to the next line. Thus, if you run
   `run.sh` again after a stop, it continues from the first case that does
   not exist.
2. If the mesh exists, it compares `maxCellSize` and the `fins` level of the
   mesh with the options of the line. It stops if they are different.
3. If the mesh directory exists but has no mesh, it stops.
4. If the mesh does not exist, it runs `mesh/Allmesh` with the options of
   the line. It stops if `Allmesh` stops with an error. The next lines that
   name this mesh use it again.
5. It runs `newCase.sh`. It stops if `newCase.sh` stops with an error.
6. It runs `./Allrun` with the default number of processes. If `Allrun`
   stops with an error, it shows a message and goes to the next line.

The cases run one after the other.

## 7. Refine the mesh near the fins (optional)

### 7.1 Command

```bash
./Allrefine          # same as ./Allrefine tip
./Allrefine tip      # cell set finTip: a band around the fin tip
./Allrefine fins     # cell set finBox: the full fin, from root to tip
```

Do this procedure in the case directory, before `./Allrun`. Alternatively,
give `--refine` to `newCase.sh`.

### 7.2 Steps that Allrefine does

1. It stops if `0/U` exists. `refineMesh` maps the fields, so a refinement
   after a run changes the results and not a clean mesh.
2. It stops if `constant/meshInfo` does not have `fins true` and `finTipR`.
3. It runs `topoSet`. `system/topoSetDict` makes the cell sets from the fin
   geometry in `constant/meshInfo`.
4. It runs `refineMesh`. Each cell in the set becomes 8 cells.
5. It removes `0/polyMesh`, which `refineMesh` writes.
6. It runs `checkMesh` and `renumberMesh`.
7. It writes the new number of cells into `nCells` in `constant/meshInfo`.

| Set | Region | Effect on y+ |
|---|---|---|
| `finTip` | Annulus from `finTipR − tipBand` to `finTipR + margin` | None. The region does not touch the body. |
| `finBox` | Cylinder of radius `finTipR + margin` | The first cell on the body below the fins is half as thick, so y+ there is approximately half. |

> **NOTE** The cell counts in the header of `common/Allrefine` are from the
> previous structured mesh. Record the counts again on the cfMesh mesh.

## 8. Results

| Data | Location |
|---|---|
| `Cd`, `Cl`, `Cm` for each 50 time steps | `postProcessing/forceCoeffs1/<start-time>/coefficient.dat` |
| Forces and moments in N, pressure and viscous parts | `postProcessing/forces1/<start-time>/force.dat`, `moment.dat` |
| Residuals of `U`, `e`, `k`, `omega` | `postProcessing/residuals/<start-time>/solverInfo.dat` |
| y+ for each wall patch | `log.rhoCentralFoam`, and the field `yPlus` in each written time |
| Mach number, pressure coefficient, wall shear stress | Fields in each written time |
| ParaView | `paraview case.foam`. `newCase.sh` makes `case.foam`. |

`<start-time>` is `0` for a new run. A continued run writes a new directory
with its start time.

`purgeWrite 2` keeps only the last two time directories.

## 9. Typical studies

### 9.1 Mach sweep

1. Write one line for each Mach number in `sweep.txt`.
2. Run `./run.sh`.
3. Measure y+ in each case. The mesh does not set y+, and the velocity
   changes by a factor of 6 from M 0.3 to M 1.8.

### 9.2 Mesh convergence

1. Change only `maxCellSize`, with a constant ratio. Do not change the
   levels. All cell sizes then change by the same ratio.
2. Use the ratio √2. The method of Celik et al. (2008) requires a ratio
   larger than 1.3. Each step multiplies the number of cells by up to
   2√2 ≈ 2.8.
3. Use the same Mach number for all meshes.

```
m03-h1  --mesh h1 --maxCellSize 1.6  --Minf 0.3
m03-h2  --mesh h2 --maxCellSize 1.13 --Minf 0.3
m03-h3  --mesh h3 --maxCellSize 0.8  --Minf 0.3
```

The first layer thickness changes with the surface cell size, so y+ changes
between the meshes. Measure y+ on each mesh. VALIDATION.md §4 gives the
analysis.

### 9.3 Angle-of-attack sweep

```bash
for a in 0 2 4 6 8; do
    ./newCase.sh ~/runs/a$a --mesh ~/meshes/base --Minf 0.8 --alpha $a --np 8
    ( cd ~/runs/a$a && ./Allrun 8 )
done
```

### 9.4 Finer fins

Use one of these two methods:

- Make a new mesh with `--finLevel 10`.
- Run `./Allrefine tip` in the case, before `./Allrun`.

## 10. Rules

- Do not put meshes in the repository. `mesh/Allmesh` makes them again.
  Put the STL and `meshDict` in the repository.
- When a change to `meshDict` changes the mesh, commit it with the number
  of cells and the `checkMesh` output of the new mesh.
- Put case directories outside the repository, outside OneDrive and outside
  `/mnt/c`.

## 11. Quick reference

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc

mesh/Allmesh ~/meshes/base                                   # mesh
./newCase.sh ~/runs/NAME --mesh ~/meshes/base --Minf 0.8 --np 8
cd ~/runs/NAME && ./Allrun 8                                 # run

./run.sh sweep.txt ~/runs ~/meshes                           # sweep

./Allrefine tip                                              # before Allrun
./Allclean                                                   # removes mesh and results
```
