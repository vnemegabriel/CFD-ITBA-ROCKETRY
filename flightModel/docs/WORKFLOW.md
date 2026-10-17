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
| `mesh/makeRocket.py` | Script that writes the geometry: the STL and its fin-edge file. MESH.md §3.3. |
| `mesh/stl/Aconcagua_biconvex117.stl` | Geometry of the rocket: the output of `makeRocket.py` with its defaults. |
| `mesh/stl/Aconcagua_biconvex117_finEdges.obj` | Leading and trailing edges of the fins, for `edgeMeshRefinement`. |
| `mesh/stl/Aconcagua.stl` | Previous geometry: 6 mm fins with bevelled edges. |
| `mesh/system/` | `meshDict`, `createPatchDict`, and the `controlDict`, `fvSchemes`, `fvSolution` that the mesh applications need. |
| `mesh/meshInfo` | Reference values and fin geometry of the mesh. |
| `case-central/` | Case template for `rhoCentralFoam`, M > 1. |
| `case-pimple/` | Case template for `rhoPimpleFoam`, M < 1. It holds only the files that differ from `case-central/`. SOLVERS.md §16. |
| `common/Allrun`, `common/Allrefine` | Scripts that `newCase.sh` copies into each case. |
| `newCase.sh` | Script that makes one case from the template and one mesh. |
| `runMesh.sh`, `meshes.txt` | Script and input file that build the meshes. |
| `runCase.sh`, `sweepSupersonic.txt`, `sweepSubsonic.txt` | Script and input files for a sweep of cases on existing meshes. |

## 2. Requirements

| Item | Requirement |
|---|---|
| OpenFOAM | ESI v2412, with cfMesh (`cartesianMesh`). |
| `jinja2` | `jinja2-cli`. Install it with `pipx install jinja2-cli`. |
| Shell tools | `bash`, `rsync`, `awk`, GNU `sed`. |
| Python 3 | Only for `mesh/makeRocket.py`. Standard library only. |
| Memory | 32 GB. `checkMesh` used 12.7 GB on a mesh of 14.2 M cells; `h3` has 10.0 M. |

Do this step before you use a script:

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
```

All scripts stop with the message `source the OpenFOAM bashrc first` if you
do not do this step.

> **CAUTION** Put the mesh directories and the case directories outside
> `/mnt/c`. That file system makes the OpenFOAM input and output several
> times slower.

## 3. Make a mesh

### 3.1 Command

```bash
mesh/Allmesh <mesh-dir> [stl-file] [--maxCellSize <m>]
             [--nLayers <n>] [--thicknessRatio <r>] [--full]
```

| Argument | Default | Function |
|---|---|---|
| `<mesh-dir>` | none, mandatory | Directory for the new mesh. It must not exist. |
| `[stl-file]` | `mesh/stl/Aconcagua_biconvex117_roll45.stl`, or `Aconcagua_biconvex117.stl` with `--full` | Geometry. MESH.md §3 gives the requirements. `<stl-file>` without `.stl` plus `_finEdges.obj` and `_finRoots.obj` must exist next to it. |
| `--maxCellSize <m>` | value in `meshDict` (1.6 m) | Sets `maxCellSize` in the copy of `meshDict`. |
| `--nLayers <n>` | value in `meshDict` (10) | Sets `nLayers` of `cone`, `walls` and `tail`. |
| `--thicknessRatio <r>` | value in `meshDict` (1.2) | Sets `thicknessRatio` of `cone`, `walls` and `tail`. |
| `--full` | half mesh | Mesh of the whole rocket, fins at 0, 90, 180 and 270°. Without it the mesh is the half z ≥ 0, cut from the full mesh with the fins at 45°, with the patch `symm` on z = 0 (MESH.md §5.1). |

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
| `<mesh-dir>/meshInfo` | A copy of `mesh/meshInfo`. `Allmesh` adds `nCells`, `maxCellSize`, `finLevel`, `nLayers`, `thicknessRatio` and `meshDir`. |
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
| `--socket <S>` | Pins the parallel run to socket `S`, one process per physical core. Refer to MEASUREMENTS.md §2. |
| `--refine tip` or `--refine fins` | Runs `./Allrefine` on the new case. Refer to §7. |
| `--init <run-dir>` | Starts from the latest time of another run instead of the uniform free stream. Refer to §5.6. |
| `--solver central` or `--solver pimple` | Template: `case-central/` (default) or `case-pimple/`. SOLVERS.md §1. |

The flow options change entries in `system/flowConditions`. The options that
you do not give keep the template values.

### 4.2 Steps that newCase.sh does

1. It stops if `--mesh` is missing, if the mesh directory has no
   `constant/polyMesh`, or if `<run-dir>` exists.
2. It copies `case-central/` into `<run-dir>`. It does not copy the mesh,
   the logs, `0/`, `processor*`, `postProcessing/` or the time directories.
   With `--solver pimple` it leaves out the files listed in
   `case-pimple/central.exclude` and then copies `case-pimple/` on top.
3. It links `<run-dir>/constant/polyMesh` to `<mesh-dir>/constant/polyMesh/`. With `--refine` it copies it. Refer to MEASUREMENTS.md §3.
4. It copies `<mesh-dir>/meshInfo` to `<run-dir>/constant/meshInfo`. This
   file replaces the `meshInfo` of the template.
5. It copies `common/Allrun` and `common/Allrefine` into `<run-dir>`.
6. It renders each `*.j2` file with the values in `config.json`, and then
   removes the `*.j2` file. Refer to SOLVERS.md §5.
7. It writes each flow option into `system/flowConditions`. It stops if the
   entry does not exist in the file.
8. It writes `--np` into `system/decomposeParDict`.
9. It runs `./Allrefine` if you give `--refine`.
10. With `--init`, it writes the path of that run into `initFrom` and into
    `run.info`.

> **NOTE** In `newCase.sh`, `--mesh` is a path. Only `runCase.sh` accepts a
> mesh name without a path.

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
   If `initFrom` exists, it maps the latest time of that run onto `0/`.
   Refer to §5.6.
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

### 5.6 Start from another run

`newCase.sh --init <run-dir>` starts the case from the solution of another
run instead of the uniform free stream. The sweep uses it from a coarse mesh
to the next finer one at the same Mach number (h1 → h2 → h3), as in the
OpenFOAM ToolChain (report §5.4): the fine mesh starts with the shock, the
wake and the boundary layer already in place, and skips the impulsive start.
The two runs can use different solvers. `Allrun` then does, after `restore0Dir`:

1. It reads the value of each `fixedValue` patch in `0/`. These are the
   values of this run: in a supersonic case, `p` on `inlet`.
2. `mapFields <run-dir> -consistent -sourceTime latestTime`. It maps the
   internal field and the patch values. The other run can have another
   mesh of the same domain and patches.
3. `changeDictionary -dict system/changeDictionaryDict.init` writes back the
   values of step 1.

The free-stream entries of the patches (`freestreamValue`, `inletValue`,
`fieldInf`) are not mapped, so they keep the values of this run. The other
run must have a time directory with the solution: in parallel, `Allrun` runs
`reconstructPar -latestTime`.

## 6. Meshes and sweeps

Meshes and cases are two steps. `runMesh.sh` builds the meshes;
`runCase.sh` makes and runs the cases on meshes that already exist.

### 6.1 Build the meshes

`meshes.txt` has one mesh on each line: a name, then the options of
`mesh/Allmesh`.

```
h1  --maxCellSize 1.6    --nLayers 10 --thicknessRatio 1.2
h2  --maxCellSize 1.1314 --nLayers 8  --thicknessRatio 1.229
```

```bash
./runMesh.sh [meshes-file] [meshes-dir]      # default: meshes.txt, ~/meshes
```

For each line:

1. If the mesh exists, it compares `maxCellSize`, the `fins` level,
   `nLayers` and `thicknessRatio` with the line. It keeps the mesh if they
   are equal and stops if they are different. It does not compare the
   geometry: after a change of the STL, remove the old meshes.
2. If the mesh directory exists but has no mesh, it stops.
3. Otherwise it runs `mesh/Allmesh` with the options of the line, and stops
   if `Allmesh` stops with an error.

### 6.2 Run the cases

A sweep file has one case on each line: a case name, then the options of
`newCase.sh`.

```
m18-h1  --mesh h1 --solver central --Minf 1.8 --pInf 78354 --Tinf 295.82
m18-h2  --mesh h2 --solver central --Minf 1.8 --pInf 78354 --Tinf 295.82  --init m18-h1
```

- Lines that start with `#`, and empty lines, have no effect.
- `--mesh` is mandatory on each line. A name without `/` is a directory in
  `<meshes-dir>`.
- `--init` is optional. A name without `/` is a directory in `<runs-dir>`.
  The run must come earlier in the file, because the cases run in order:
  `m12-h1 --mesh h1 --Minf 1.2 ... --init m18-h1`.

```bash
./runCase.sh <sweep-file> [runs-dir] [meshes-dir]   # default: ~/runs, ~/meshes
```

For each line:

1. If `<runs-dir>/<name>` exists, it goes to the next line. Thus, if you run
   `runCase.sh` again after a stop, it continues from the first case that
   does not exist.
2. If the mesh does not exist, it stops. Build it with `runMesh.sh`.
3. It runs `newCase.sh`. It stops if `newCase.sh` stops with an error.
4. It runs `./Allrun` with the default number of processes. If `Allrun`
   stops with an error, it shows a message and goes to the next line.

The cases run one after the other. To use both sockets, run two sweep files
at the same time, one with `--socket 0` and one with `--socket 1`
(MEASUREMENTS.md §2). Each chain of `--init` must stay in one file.

A sweep runs for days. Start it with `nohup` and `&`:

```bash
nohup ./runCase.sh sweepSupersonic.txt > ~/log.sweepSupersonic 2>&1 &
tail -f ~/log.sweepSupersonic
```

- `nohup` keeps `runCase.sh` alive when the terminal closes or the ssh
  session drops. Without it, the shell sends SIGHUP and the sweep stops.
- `> <log> 2>&1` writes the output and the errors into the log. Without it,
  `nohup` writes `nohup.out` in the current directory.
- `&` runs it in the background, so you can start the second sweep from the
  same terminal.

`tmux` or `screen` do the same and let you attach to the session again.

To stop a sweep, kill `runCase.sh` (`pkill -f 'runCase.sh sweepSupersonic'`).
The solver of the current case goes on until it ends; no new case starts.
To stop the solver too, kill its `mpirun`.

> **CAUTION** `runCase.sh` reads the sweep file line by line while it runs.
> Do not edit the file in place until it ends: the next line is read from the
> old byte offset of the new content. Stop `runCase.sh` (the solver can go
> on), edit, and run it again; it skips the cases that exist.

> **CAUTION** If your home directory is in `/mnt/c`, give a `runs-dir` and a
> `meshes-dir` outside it.

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

## 8. Results

| Data | Location |
|---|---|
| `Cd`, `Cl`, `Cm` for each 50 time steps | `postProcessing/forceCoeffs1/<start-time>/coefficient.dat` |
| Forces and moments in N, pressure and viscous parts | `postProcessing/forces1/<start-time>/force.dat`, `moment.dat` |
| Residuals of `U`, `e`, `k`, `omega` | `postProcessing/residuals/<start-time>/solverInfo.dat` |
| y+ for each wall patch | `postProcessing/yPlus/<start-time>/yPlus.dat`, and the field `yPlus` in the last written time |
| Mach number, pressure coefficient, wall shear stress | Fields in each written time |
| p, U, T, ρ, Mach, Cp, wall shear stress, y+ on the symmetry plane and on the rocket skin | `postProcessing/surfaces/<end-time>/{symm,skin}.vtp`, written at the end of the run |
| ParaView | `paraview case.foam`. `newCase.sh` makes `case.foam`. |

`<start-time>` is `0` for a new run. A continued run writes a new directory
with its start time.

`purgeWrite 1` keeps only the last time directory. MEASUREMENTS.md lists all the outputs, `run.info` and the notification.

### 8.1 Monitor a run

```bash
./post/monitor.sh <run-dir> [seconds]     # default 30 s
```

It plots `Cd` and `Cl` against the iteration and the initial residuals in log
scale, and draws them again every `[seconds]`. It joins the directories of
a continued run. Without a display, for example over ssh, it writes
`<run-dir>/monitor.png`. It needs gnuplot (`sudo apt install gnuplot`).

### 8.2 Table of a sweep

```bash
./post/sweepTable.sh [runs-dir] > sweep.csv    # default ~/runs
```

One line for each run. A run that has not started has empty fields.

| Columns | Source |
|---|---|
| `case`, `mesh`, `nCells` | `run.info` |
| `solver` | `application` in `system/controlDict` |
| `Minf`, `pInf`, `Tinf`, `alpha` | `system/flowConditions` |
| `status`, `steps`, `tEnd`, `wallSeconds` | Last record of `run.info` |
| `Cd`, `Cl`, `CmPitch` | Last line of `forceCoeffs1` |
| `window`, `CdMean`, `ClMean`, `CmPitchMean` | Mean over the last `coeffs_window` iterations (`config.json`), the window of the stop criterion |
| `resUx` … `resOmega` | Last initial residuals of `residuals`; `resP` is empty for `rhoCentralFoam`, which solves no pressure equation |
| `yPlusMin`, `yPlusMax`, `yPlusAvgWalls`, `yPlusAvgFins` | `yPlus` at the last time: over all the wall patches, and the mean on `walls` and on `fins` |

With a continued run, the coefficients and residuals come from the last
start-time directory only.

## 9. Typical studies

### 9.1 Mach sweep

1. Write one line for each Mach number in a sweep file.
2. Run `./runCase.sh <sweep-file>`.
3. Measure y+ in each case. The mesh does not set y+, and the velocity
   changes by a factor of 6 from M 0.3 to M 1.8.

### 9.2 Mesh convergence

1. Change only `maxCellSize`, with a constant ratio. Do not change the
   levels. All cell sizes then change by the same ratio.
2. Use the ratio √2. The method of Celik et al. (2008) requires a ratio
   larger than 1.3. Each step multiplies the number of cells by up to
   2√2 ≈ 2.8. The mesh does not refine uniformly: measured, N grows ×1.95
   from `h1` to `h2` and ×2.34 from `h2` to `h3`, so r = (N_fine/N_coarse)^(1/3)
   is 1.25 and 1.33. Use these values, not √2, in VALIDATION.md §4.2.
3. Use the same Mach number for all meshes.

```
# meshes.txt
h0  --maxCellSize 2.2627 --nLayers 12 --thicknessRatio 1.1878
h1  --maxCellSize 1.6    --nLayers 10 --thicknessRatio 1.2
h2  --maxCellSize 1.1314 --nLayers 8  --thicknessRatio 1.229
h3  --maxCellSize 0.8    --nLayers 7  --thicknessRatio 1.2

# sweep file
m03-h0  --mesh h0 --solver pimple --Minf 0.3
m03-h1  --mesh h1 --solver pimple --Minf 0.3  --init m03-h0
m03-h2  --mesh h2 --solver pimple --Minf 0.3  --init m03-h1
m03-h3  --mesh h3 --solver pimple --Minf 0.3  --init m03-h2
```

4. Keep the first cell the same on every mesh. `nLayers` and
   `thicknessRatio` hold the layer rule of MESH.md §4.4 at 0.241 mm, which
   `cartesianMesh` makes about 0.10 mm. Measure it on each new mesh:

   ```bash
   ./post/firstLayer.py ~/meshes/h2
   ```

   The median t₁ on `walls` must agree with that of `h1` (0.100 mm). If not,
   y+ changes between the meshes and the study mixes the wall treatment with
   the discretization error.
5. `h0` is not part of the study. It is a coarser mesh to start the `h1`
   runs (`--init`, §5.6); each run then starts from the coarser one.

MESH.md §4.4 gives the y+ at each Mach number. Measure y+ on each run to
confirm. VALIDATION.md §4 gives the analysis.

### 9.3 Angle-of-attack sweep

```bash
for a in 0 2 4 6 8; do
    ./newCase.sh ~/runs/a$a --mesh ~/meshes/base --Minf 0.8 --alpha $a --np 8
    ( cd ~/runs/a$a && ./Allrun 8 )
done
```

### 9.4 Finer fins

Use one of these two methods:

- Raise the levels of `fins`, `finEdges` and `finRoots` in `mesh/system/meshDict` and make a new mesh.
- Run `./Allrefine tip` in the case, before `./Allrun`.

### 9.5 The sweep of record: 3 meshes × 5 Mach numbers × 3 angles of attack

`meshes.txt` holds the three meshes of the study plus `h0` to start them,
and two sweep files the 66 runs: 5 Mach numbers × α = 0, 4, 8° × 4 meshes,
plus 6 `pimple` start runs. `sweepSupersonic.txt` has M 1.2 and 1.8 (30
runs) on socket 0, `sweepSubsonic.txt` has M 0.3, 0.6 and 0.8 (36 runs) on
socket 1, each with 14 processes, one per physical core. Each Mach number and angle is a chain
h0 → h1 → h2 → h3, each run started from the previous one (`--init`). M 0.3,
0.6 and 0.8 run `pimple` on all the meshes. M 1.2 and 1.8 start with a
`pimple` run on `h0` (`m18-a04-h0-pimple`, ...) and run `central` from it
on, because `central` from the uniform free stream crashes at the base
(SOLVERS.md §17). The runs are in order of mesh, so all the `h0` results
come first.

`--alpha` turns the free stream in the x–y plane, the plane of the
half-mesh symmetry; the mesh is the same for all the angles. The wake
regions of `meshDict` (`baseWake`, `nearWake`, `midWake`) are cones on the x
axis, but the wake follows the free stream and moves tan α per metre
towards +y:

| Region | Radius | Length behind the base | Offset at 4° | Offset at 8° |
|---|---|---|---|---|
| `baseWake` | 0.08 m | 0.36 m | 0.03 m | 0.05 m |
| `nearWake` | 0.14 m | 0.85 m | 0.06 m | 0.12 m |
| `midWake` | 0.26 m | 2.1 m | 0.15 m | 0.30 m |

At 4° the wake stays inside the three regions. At 8° it leaves `midWake`
about 1.85 m behind the base. The sweep had α = 16° until 2026-10-17; it was
replaced by 4° because the TeleMega puts α in the coast below 0.2°
(`flightData/`), and at 16° the wake left the regions 0.4 m behind the base.

To make the same results on another
computer:

1. Do the steps of §2.
2. Optional: write the geometry again and compare it with the file in the
   repository. The two files must be equal.

   ```bash
   python3 mesh/makeRocket.py /tmp/check.stl && cmp /tmp/check.stl mesh/stl/Aconcagua_biconvex117.stl
   ```

3. Make the meshes and examine them:

   ```bash
   ./runMesh.sh
   ```

   The meshes must agree with this table:

   | Mesh | `maxCellSize` | Cells | Negative volumes | Max skewness | Max non-orthogonality | t₁ on `walls` |
   |---|---|---|---|---|---|---|
   | `h0` | 2.2627 | 1 074 798 | 0 | 2.37 | 65.7 | 0.097 mm |
   | `h1` | 1.6 | 2 198 939 | 0 | 1.72 | 66.0 | 0.100 mm |
   | `h2` | 1.1314 | 4 285 770 | 0 | 1.58 | 65.6 | 0.096 mm |
   | `h3` | 0.8 | 10 020 768 | 0 | 1.41 | 65.9 | 0.101 mm |

4. Run the two files at the same time:

   ```bash
   nohup ./runCase.sh sweepSupersonic.txt > ~/log.sweepSupersonic 2>&1 &
   nohup ./runCase.sh sweepSubsonic.txt   > ~/log.sweepSubsonic   2>&1 &
   ```

   On a computer with one socket of 8 cores (nemepad, while newton is out
   of service) use `sweepSubsonic.nemepad.txt` and
   `sweepSupersonic.nemepad.txt`: the same runs with `--np 8` and no
   `--socket`. Run one file, then the other, never both at once. On `h0`,
   `m08-a00-h0` takes ~2.7 s per step with 8 processes (clean intervals of the
   sweep runs, 2026-10-17), against 1.77 s with 14 processes on one socket of
   newton: about 1.5 times slower per step. The memory is about 3 GB per
   million cells: `h3` does not fit in the 28 GB of nemepad.

The flight conditions are those of the coast phase (t > 4 s) of the IREC
flight in Midland, with the atmosphere of that day, not the ISA. The data is
in `flightData/AconcaguaIREC_TeleMega.csv`, and `flightData/cdFlight.py`
gets everything from it:

- The pressure at the pad (896 m MSL) is the mean before launch, 90565 Pa.
- The temperature profile is fitted to the GPS height above the pad against
  the measured pressure, on the descent: 305.0 K at the pad and 7.2 K/km,
  about ISA + 20 K.
- The speed is the axial accelerometer integrated from t = 18 s, where the
  baro is clean. The Mach number comes from that speed and the fitted
  temperature; `pInf` and `Tinf` come from the integrated height.

| M | t [s] | Altitude MSL [m] | `pInf` [Pa] | `Tinf` [K] |
|---|---|---|---|---|
| 1.8 | 4.1 | 2171 | 78354 | 295.82 |
| 1.2 | 8.5 | 4347 | 60547 | 280.15 |
| 0.8 | 14.6 | 6289 | 47509 | 266.17 |
| 0.6 | 19.8 | 7437 | 40911 | 257.90 |
| 0.3 | 29.0 | 8679 | 34613 | 248.96 |

Until 2026-10-17 they were the ISA values at the altitude of the TeleMega
(e.g. M 0.8: 49800 Pa, 251.72 K); the runs in `~/runs.sin-medicion` use those.

## 10. Rules

- Do not put meshes in the repository. `mesh/Allmesh` makes them again.
  Put the STL and `meshDict` in the repository.
- When a change to `meshDict` changes the mesh, commit it with the number
  of cells and the `checkMesh` output of the new mesh.
- Put case directories outside the repository and outside `/mnt/c`.

## 11. Quick reference

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc

mesh/Allmesh ~/meshes/base                                   # mesh
./newCase.sh ~/runs/NAME --mesh ~/meshes/base --Minf 0.8 --np 8
cd ~/runs/NAME && ./Allrun 8                                 # run

./runMesh.sh meshes.txt ~/meshes                             # meshes
./runCase.sh sweepSupersonic.txt ~/runs ~/meshes             # sweep

./Allrefine tip                                              # before Allrun
./Allclean                                                   # removes mesh and results
```
