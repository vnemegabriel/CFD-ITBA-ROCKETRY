# CFD-ITBA-ROCKETRY

Aerodinámica del cohete Aconcagua en OpenFOAM v2412. Malla con cfMesh, casos
para subsónico, transónico y supersónico, y un barrido de Mach en serie.

```
mesh/         Aconcagua.stl, Allmesh, meshDict — la malla
common/       Allrun, Allrefine — los scripts que va a usar cada corrida
case-*/       las tres plantillas de caso
newCase.sh    arma una corrida: plantilla + malla + condiciones de vuelo
run.sh        corre todo sweep.txt en serie
sweep.txt     un punto de Mach por línea
docs/         cómo se usa y por qué
```

## Arranque

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc

cd mesh && ./Allmesh && cd ..                 # ~5 min, 862 k celdas
./newCase.sh ~/runs/prueba --regime sub --Uinf 68 --np 8
cd ~/runs/prueba && ./Allrun 8
```

El barrido entero:

```bash
./run.sh
```

Los directorios de corrida van **fuera de OneDrive y fuera de `/mnt/c`**: las
dos cosas hacen la E/S de OpenFOAM varias veces más lenta.

## Qué leer

| | |
|---|---|
| [WORKFLOW.md](docs/WORKFLOW.md) | el flujo completo, de la geometría al Cd |
| [SOLVERS.md](docs/SOLVERS.md) | qué solver usa cada régimen y por qué |
| [VALIDATION.md](docs/VALIDATION.md) | contra qué comparar los resultados |
| [TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md) | lo que ya falló y cómo se arregló |

## Requisitos

OpenFOAM v2412 de ESI, que ya trae cfMesh (`cartesianMesh`). Bajo Windows,
en WSL. Nada de Python.

## Régimen

| régimen | Mach | solver |
|---|---|---|
| `sub` | < 0.3 | `simpleFoam` |
| `trans` | 0.3 – 1.2 | `rhoSimpleFoam` |
| `super` | > 1.2 | `rhoCentralFoam` |

## Estado

La malla actual: 862 k celdas, 95 % hexaedros, sin volúmenes negativos,
no-ortogonalidad máxima 78.6 y media 5.7.

Pendiente: medir y+ en cada punto del barrido y ajustar `nLayers`. Las aletas
van con 3 capas porque el borde de ataque es un filo y más capas enredan la
extrusión — está medido en TROUBLESHOOTING.
