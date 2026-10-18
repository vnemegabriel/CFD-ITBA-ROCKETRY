# CFD-ITBA-ROCKETRY

Aerodynamics of the Aconcagua rocket with OpenFOAM v2412: cfMesh meshes,
`rhoPimpleFoam` for M < 1 and `rhoCentralFoam` for M > 1, Mach sweeps in
series, and a comparison with the flight data of the IREC flight in Midland.

**User guide:** [`flightModel/docs/guide/aconcagua-user-guide.pdf`](flightModel/docs/guide/aconcagua-user-guide.pdf)
(source in `flightModel/docs/guide/`). It is the only user documentation.

```
flightModel/mesh/         makeRocket.py, stl/, Allmesh, meshInfo, system/   the mesh
flightModel/common/       Allrun, Allrefine, machineLog.py                  copied into each run
flightModel/case-central/ rhoCentralFoam template, M > 1
flightModel/case-pimple/  rhoPimpleFoam template, M < 1 (only what differs)
flightModel/newCase.sh    one run: template + mesh + flight condition
flightModel/runMesh.sh    the meshes of meshes.txt
flightModel/runCase.sh    one sweep file in series, on meshes that exist
flightModel/post/         monitor.sh, sweepTable.sh, firstLayer.py
flightModel/flightData/   TeleMega data and cdFlight.py (flight Cd)
flightModel/docs/guide/   the user guide (LaTeX)
```

## First run

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
cd flightModel

mesh/Allmesh ~/meshes/h1
./newCase.sh ~/runs/test --mesh ~/meshes/h1 --solver pimple --Minf 0.8 --np 8
cd ~/runs/test && ./Allrun 8
```

Put the mesh and run directories outside `/mnt/c`. The user guide gives the
requirements, the preparation of the computer for a sweep, and the procedures.

`fin_flutter/` is a separate project.
