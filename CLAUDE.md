# CLAUDE.md — Aconcagua CFD

Aerodynamics of the Aconcagua rocket with OpenFOAM v2412. cfMesh meshes, a
`rhoPimpleFoam` template for M < 1 and a `rhoCentralFoam` template for M > 1
that starts from a pimple run at the same Mach number, Mach sweeps in series.

Reply in the language of the person who asks; the maintainers use Rioplatense
Spanish. CFD jargon stays in English (mesh, solver, y+, wall function,
skewness). Code, comments, scripts and documentation are in English. The
documentation is the LaTeX user guide in `flightModel/docs/guide/`, in
Simplified Technical English (ASD-STE100) with the structure of ISO/IEC/IEEE
26514. A change of a procedure changes the guide in the same commit: compile it
with `make` there (lualatex) and commit the PDF too. The guide gives the method
only; measured results go to the technical report in `flightModel/docs/report/`.

**Few comments in the code.** The documentation explains. A comment is
justified only when the code cannot say it — for example, that the six numbers
of `surfaceGenerateBoundingBox` are metres and not multiples of the bounding
box.

---

## The objective

**Cd against Mach number: 0.3, 0.6, 0.8, 1.2, 1.8.** What the mesh must
resolve, in this order:

1. The shocks, for M ≥ 0.8.
2. The base wake, a large part of the Cd of a rocket.
3. The wall friction, which depends on y+.

The sweep runs at α = 0, 4 and 8°, set with the velocity field (`--alpha`): the
mesh is the same. The wake cones of `meshDict` follow the x axis; at 8° the
wake leaves `midWake` about 1.85 m behind the base. Do not spend cells on the
fin-tip vortex.

**One set of meshes for the whole sweep:** h1, h2, h3 for the convergence
study (only `maxCellSize` changes, ratio √2) and h0, coarser, which always
starts the h1 runs (`newCase.sh --init`). Half mesh, with the symmetry plane at
z = 0 and the fins at 45°: `Allmesh` meshes the whole rocket and cuts it at
z = 0 with `subsetMesh`. The layers keep the first body cell equal on all
meshes; `nutUSpaldingWallFunction` covers the buffer layer. Measure t1 on each
new mesh with `post/firstLayer.py`. Do not make meshes per Mach number. Measure
the real y+ on each run.

---

## The structure

```
flightModel/mesh/         makeRocket.py, stl/, Allmesh, meshInfo, system/
flightModel/common/       Allrun, Allrefine, machineLog.py
flightModel/case-central/ rhoCentralFoam template, M > 1
flightModel/case-pimple/  rhoPimpleFoam template, M < 1 and start of M > 1: only what differs from central
flightModel/newCase.sh    template + mesh + conditions -> one run
flightModel/runMesh.sh    meshes.txt: builds the meshes that are missing
flightModel/runCase.sh    one sweep*.txt in series, on existing meshes; --np and --socket
                          are its options, not lines of the sweep file
flightModel/post/         monitor.sh (live gnuplot), sweepTable.sh (CSV),
                          firstLayer.py (t1 and y of the first cell)
flightModel/flightData/   TeleMega telemetry (dot decimal), Aconcagua.ork, cdFlight.py (flight Cd)
flightModel/docs/guide/   the user guide (method only): LaTeX, Mermaid diagrams rendered with Pele
flightModel/docs/report/  the technical report (results), for the thesis
```

The Python is `mesh/makeRocket.py`, which writes the STL, `mesh/snapPlane.py`,
which flattens the cut of the half mesh, `post/firstLayer.py`, which measures
the first layer, `common/machineLog.py`, which records the computer during each
run, and `flightData/cdFlight.py`. All use only the standard library, except
matplotlib in `cdFlight.py`. The mesh is defined in `mesh/system/meshDict`, and
everything there is a **level**: an integer that halves a cell, local to a
patch or a box, not coupled to anything else.

---

## Rules

- **Measure before you propose.** `mesh/Allmesh <dir>` gives cells, negative
  volumes, non-orthogonality and skewness in a few minutes. A change without a
  number next to it is not finished.
- **`checkMesh` with negative volumes = unusable mesh.** The other
  `-allGeometry` checks that fail today (small determinant in layer cells, face
  tets, concave cells) are of another category and do not block.
- **The fins use `nLayers 3`.** More tangles the extrusion. Measured on the old
  6 mm fin with a sharp leading edge; the current fin is 11.7 mm thick with
  1 mm edge radii, and it was not measured again.
- **The fin edges have a radius.** With sharp edges, the layers of the two faces
  collide at the edge. Guide, chapter Geometry.
- **The fin root has a fillet** (`makeRocket.py --finRootRadius`). Without it the
  layers of the fin and of the body close at the corner. Guide, chapter Geometry.
- Before you cut or change something, grep `docs/guide/chapters/` as well as the
  code.
- One task per session. One commit per task.

## Known traps

- **With `localEuler` the Courant number in the log means nothing.**
  rhoCentralFoam calculates it with `deltaT 1` of `controlDict` (it gives ~1e7)
  and rhoPimpleFoam does not write it. Read `Flow time scale min/max`. This
  misreading kept the central template on a global `Euler` step until
  2026-10-12. Guide, "Local time step".
- **Scripts that calculate with `awk` set `LC_ALL=C`.** In a locale with a
  decimal comma (`es_ES`), `awk` reads and writes numbers with a comma:
  `Allmesh` once gave `11,775` to `surfaceGenerateBoundingBox` and made a wrong
  domain.
- **Do not put Python scripts in bash here-documents**: a `\n` inside an f-string
  becomes a real new line and breaks the file.
- **`runApplication` returns 0 even when the application dies** with FOAM FATAL
  ERROR, so `set -e` does not catch it. `mesh/Allmesh` and `Allrun` check each
  log for this reason.
- **`mapFields` leaves the cells it cannot find with the values of `0.orig`.**
  `0.orig` has no `rho`, so those cells got rho = p∞ and the h1 runs at α ≠ 0
  died at step 2. `Allrun` deletes `0/rho` after the mapping.
- A `*/` inside the banner comment of an OpenFOAM dictionary closes the comment
  early, and the parser reads garbage.
- Do not edit a shell script while it runs: `sh` reads it as it executes it.
