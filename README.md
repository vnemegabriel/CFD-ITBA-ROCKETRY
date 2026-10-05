# CFD-ITBA-ROCKETRY

Aerodinámica del cohete Aconcagua en OpenFOAM v2412. Malla con cfMesh, un solo
solver (`rhoCentralFoam`) para todo el barrido de Mach, en serie.

```
mesh/         stl/Aconcagua.stl, Allmesh, meshDict — la malla
common/       Allrun, Allrefine — los scripts que va a usar cada corrida
case-central/ la plantilla de caso
newCase.sh    arma una corrida: plantilla + malla + condiciones de vuelo
run.sh        corre todo sweep.txt en serie
sweep.example.txt  plantilla de barrido; se copia a sweep.txt, que no se versiona
docs/         cómo se usa y por qué
```

## Arranque

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc

mesh/Allmesh ~/meshes/base                    # ~2 min, 2.1 M celdas
./newCase.sh ~/runs/prueba --mesh ~/meshes/base --Minf 0.8 --np 8
cd ~/runs/prueba && ./Allrun 8
```

El barrido entero:

```bash
./run.sh
```

Las corridas van en `~/runs` y las mallas en `~/meshes`, **fuera de la repo**:
pesan gigas y no se versionan.

## Qué leer

| | |
|---|---|
| [WORKFLOW.md](docs/WORKFLOW.md) | el flujo completo, de la geometría al Cd |
| [SOLVERS.md](docs/SOLVERS.md) | el solver, los modelos y por qué |
| [VALIDATION.md](docs/VALIDATION.md) | contra qué comparar los resultados |
| [TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md) | lo que ya falló y cómo se arregló |

## Requisitos

OpenFOAM v2412 de ESI, que ya trae cfMesh (`cartesianMesh`). Bajo Windows,
en WSL. La plantilla necesita `jinja2-cli`
(`pip install jinja2-cli`) para renderizar sus `*.j2`; ver
[SOLVERS.md §2.3](docs/SOLVERS.md).


