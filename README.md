# CFD-ITBA-ROCKETRY

Aerodynamics of the Aconcagua rocket with OpenFOAM v2412: cfMesh meshes,
`rhoPimpleFoam` for M < 1 and `rhoCentralFoam` for M > 1, Mach sweeps in
series, and a comparison with the data of the IREC flight in Midland.

| Document | Contents |
|---|---|
| [User guide](flightModel/docs/guide/aconcagua-user-guide.pdf) | The method: installation, geometry, mesh, cases, sweeps, measurement, validation. Start here. |
| [Technical report](flightModel/docs/report/aconcagua-cfd-report.pdf) | The results, and the reference values to check a replication. |
| [CLAUDE.md](CLAUDE.md) | Rules for people and agents who change the repository. |

```
container/Dockerfile      OpenFOAM v2412 + cfMesh + tools, for any computer
flightModel/mesh/         makeRocket.py, stl/, Allmesh, meshInfo, system/   the mesh
flightModel/common/       Allrun, Allrefine, machineLog.py                  copied into each run
flightModel/case-central/ rhoCentralFoam template, M > 1
flightModel/case-pimple/  rhoPimpleFoam template, M < 1 (only what differs)
flightModel/newCase.sh    one run: template + mesh + flight condition
flightModel/runMesh.sh    the meshes of meshes.txt
flightModel/runCase.sh    one sweep file in series, on meshes that exist
flightModel/post/         monitor.sh, sweepTable.sh, firstLayer.py
flightModel/flightData/   TeleMega telemetry, OpenRocket file, cdFlight.py (flight Cd)
flightModel/docs/         the user guide and the technical report (LaTeX)
```

## Quick start (any computer, with Docker)

```bash
git clone https://github.com/vnemegabriel/CFD-ITBA-ROCKETRY.git
cd CFD-ITBA-ROCKETRY
docker build -t aconcagua-cfd -f container/Dockerfile .
docker run --rm -it -v "$PWD":/work -v aconcagua-data:/data aconcagua-cfd
```

Inside the container, make the coarsest mesh, h0, and run one case. h0 always
starts the chains of a study, so this is also the first run of the sweep:

```bash
mesh/Allmesh ~/meshes/h0 $(grep '^h0 ' meshes.txt | cut -d' ' -f2-)
./newCase.sh ~/runs/m08-a00-h0 --mesh ~/meshes/h0 --solver pimple --alpha 0 \
    --Minf 0.8 --pInf 47509 --Tinf 266.17
cd ~/runs/m08-a00-h0 && ./Allrun
```

Compare the number of cells and Cd with the reference values of the technical
report. The native installation on Ubuntu or WSL 2 is in the user guide,
chapter "Installation and setup".

`fin_flutter/` is a separate project.
