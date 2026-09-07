# CFD-ITBA-ROCKETRY

Repositorio del equipo ITBA ROCKETRY del Instituto Tecnológico de Buenos Aires.

Malla estructurada 100 % hexaédrica del cohete **Aconcagua** (punta, cuerpo,
boattail, base y aletas), generada con scripts de Python sobre gmsh, y el caso
OpenFOAM que la corre. Todo se controla desde la línea de comandos; nada se
hace a mano en una GUI.

```
CFD-ITBA-ROCKETRY/
├── mesh/        generador de malla.  python build.py  es el único comando
│   ├── meshParams.py      parámetros por defecto (el único archivo que se edita)
│   ├── aconcaguaGeom.py   geometría de registro del cohete
│   ├── presets/           smoke, coarse, medium, fine, fintip, wake_unsteady
│   └── output/            .msh + sidecars (.meshInfo, .params.py); no se versiona
├── case/        plantilla OpenFOAM v2412 (simpleFoam + kOmegaSST).  NO se corre acá
│   ├── system/flowConditions   Uinf, alpha, beta, nu ...  lo único que se edita
│   ├── Allmesh / Allrefine / Allrun / Allclean
│   └── 0.orig/ constant/ system/
├── newCase.sh   plantilla + malla + condiciones  ->  directorio de corrida
└── docs/        empezar por WORKFLOW.md
```

## Quick start

En Windows la malla se construye con el Python nativo; OpenFOAM corre en WSL
(Ubuntu, `/usr/lib/openfoam/openfoam2412`). En Linux todo corre en la misma
shell.

```bash
# 1. malla (Windows o Linux).  ~45 s, ~530 k celdas por cuadrante
cd mesh
python build.py --plan                    # ver el plan de celdas sin construir nada
python build.py --preset coarse           # -> output/aconcagua_quarter_s3.msh

# 2. caso (shell con OpenFOAM cargado)
source /usr/lib/openfoam/openfoam2412/etc/bashrc
cd ..
./newCase.sh ~/runs/prueba mesh/output/aconcagua_quarter_s3.msh
cd ~/runs/prueba && ./Allrun 8

# 3. resultados
tail postProcessing/forceCoeffs1/0/coefficient.dat
```

Ángulo de ataque: la malla de un cuarto tiene dos planos de simetría y sólo
sirve para flujo axial. Para `alpha` se construye media malla, para `alpha` y
`beta` (o balanceo) la malla completa:

```bash
python build.py --preset medium --sector half
./newCase.sh ~/runs/a05 mesh/output/aconcagua_half_s2.msh --alpha 5
```

## Qué leer

| Documento | Para qué |
|---|---|
| [docs/WORKFLOW.md](docs/WORKFLOW.md) | el flujo de trabajo completo: construir, convertir, correr, chequear. **Empezar acá.** |
| [docs/PARAMETERS.md](docs/PARAMETERS.md) | referencia de cada parámetro de la malla y del caso |
| [docs/SECTORS_AND_AOA.md](docs/SECTORS_AND_AOA.md) | cuarto / mitad / completa, ángulo de ataque, cómo se ensambla |
| [docs/GEOMETRY.md](docs/GEOMETRY.md) | la geometría de registro y cómo cambiarla (otro cohete, otras aletas) |
| [docs/MESH_DESIGN.md](docs/MESH_DESIGN.md) | por qué la topología es la que es; calidad; restricciones del toolchain (en inglés) |
| [docs/TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md) | errores conocidos y qué significan |
| [case/README.md](case/README.md) | el caso OpenFOAM archivo por archivo |

## Requisitos

- Python ≥ 3.10 con `numpy`, `gmsh` (`pip install gmsh`) y, recomendado, `scipy`.
  En Linux/WSL gmsh necesita además `sudo apt install libglu1-mesa`.
- OpenFOAM v2412 (ESI). Los scripts usan `foamDictionary`, `gmshToFoam`,
  `checkMesh`, `renumberMesh`, `potentialFoam`, `simpleFoam`.

Las mallas `.msh` no se versionan (un cuarto fino son 520 MB): se regeneran
desde su `.params.py` en segundos o minutos. Mantener `mesh/output/` y los
directorios de corrida **fuera de OneDrive** y, en WSL, fuera de `/mnt/c`.

## Pendiente

- Calculadora de fin flutter
