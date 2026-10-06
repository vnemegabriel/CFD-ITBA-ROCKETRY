# CFD-ITBA-ROCKETRY

Aerodinámica del cohete Aconcagua en OpenFOAM v2412. Malla con cfMesh, un solo
solver (`rhoCentralFoam`) para todo el barrido de Mach, en serie.

```
flightModel/mesh/         makeRocket.py, stl/, Allmesh, meshDict — la malla
flightModel/common/       Allrun, Allrefine — los scripts que va a usar cada corrida
flightModel/case-central/ la plantilla de caso
flightModel/newCase.sh    arma una corrida: plantilla + malla + condiciones de vuelo
flightModel/run.sh        corre todo sweep.txt en serie
flightModel/sweep.txt     un punto de Mach por línea
flightModel/docs/         cómo se usa y por qué
```

## Arranque

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc

flightModel/mesh/Allmesh ~/meshes/base                    # ~2 min, 2.1 M celdas
flightModel/newCase.sh ~/runs/prueba --mesh ~/meshes/base --Minf 0.8 --np 8
cd ~/runs/prueba && ./Allrun 8
```

El barrido entero:

```bash
flightModel/run.sh
```

Los directorios de corrida van **fuera de OneDrive y fuera de `/mnt/c`**: las
dos cosas hacen la E/S de OpenFOAM varias veces más lenta.

## Qué leer

| | |
|---|---|
| [WORKFLOW.md](flightModel/docs/WORKFLOW.md) | el flujo completo, de la geometría al Cd |
| [SOLVERS.md](flightModel/docs/SOLVERS.md) | el solver, los modelos y por qué |
| [VALIDATION.md](flightModel/docs/VALIDATION.md) | contra qué comparar los resultados |
| [TROUBLESHOOTING.md](flightModel/docs/TROUBLESHOOTING.md) | lo que ya falló y cómo se arregló |

## Requisitos

OpenFOAM v2412 de ESI, que ya trae cfMesh (`cartesianMesh`). Bajo Windows,
en WSL. La plantilla necesita `jinja2-cli`
(`pip install jinja2-cli`) para renderizar sus `*.j2`; ver
[SOLVERS.md §2.3](flightModel/docs/SOLVERS.md).


