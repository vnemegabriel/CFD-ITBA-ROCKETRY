# MEASUREMENTS.md — qué se mide en cada corrida y dónde queda

Cada caso que arma `newCase.sh` mide sus propios coeficientes, su costo y su
estado, y lo deja en archivos chicos. La idea es poder comparar contra otros
solvers (Nektar++) sin guardar gigas por corrida.

## 1. Qué se extrae

### Coeficientes y fuerzas, cada 50 pasos

| Archivo | Contenido |
|---|---|
| `postProcessing/forceCoeffs1/<t0>/coefficient.dat` | `Cd`, `Cl`, `CmPitch`, `CmRoll`, `CmYaw`, `Cs`, y cada uno partido en eje delantero `(f)` y trasero `(r)` |
| `postProcessing/forces1/<t0>/force.dat`, `moment.dat` | Fuerza y momento en N y N·m, total, de presión y viscoso |
| `postProcessing/residuals/<t0>/solverInfo.dat` | Residuos de `U`, `e`, `k`, `omega` |
| `postProcessing/yPlus/<t0>/yPlus.dat` | y+ mínimo, máximo y medio por patch, en cada tiempo escrito |

Referencias: `Aref = π·rBody²`, `lRef = diámetro`, `rhoInf` y `Uinf` de
`system/flowConditions`. α = 0, así que `Cl`, `Cs` y los momentos deberían
dar ≈ 0; sirven de control de simetría.

### Campos, sólo en el último tiempo escrito

`p`, `T`, `U`, `rho`, `k`, `omega`, `nut`, `alphat` (los del solver), más
`Ma`, `static(p)_coeff` (Cp), `yPlus` y `wallShearStress`.

No se escriben y salen después de `U` o `T` con `postProcess`, sin volver a
correr:

| Magnitud | Cómo |
|---|---|
| Vorticidad | `postProcess -func vorticity -latestTime` |
| Criterio Q | `postProcess -func Q -latestTime` |
| Presión y temperatura totales | function object `pressure` con `mode total`; T₀ = T·(1 + (γ−1)/2·Ma²) |
| Energía interna específica | e = cv·T, con cv = R/(γ−1) = 717.6 J/(kg K) |

## 2. Costo de la corrida

### `postProcessing/timeInfo/<t0>/timeInfo.dat`

Una línea cada 50 pasos: tiempo físico, CPU y wall acumulados, y CPU y wall
de los últimos 50 pasos. El deltaT sale de dos líneas seguidas de
`coefficient.dat`: Δt/50.

### `run.info`

`newCase.sh` lo crea con los datos del caso. `Allrun` agrega un bloque
`run_<fecha>` por cada vez que corre (una continuación agrega otro). Es un
dict de OpenFOAM, así que se lee con `foamDictionary`.

| Entrada | Quién | Qué |
|---|---|---|
| `case`, `created` | `newCase.sh` | Nombre y fecha |
| `commit` | `newCase.sh` | `git describe --dirty` de la repo; `-dirty` = había cambios sin commitear |
| `mesh`, `nCells` | `newCase.sh` | Malla enlazada y número de celdas |
| `options` | `newCase.sh` | La línea de opciones tal cual |
| `host`, `cpu`, `np`, `cpuset` | `Allrun` | Máquina, procesos MPI y CPUs permitidos (`OMPI_MCA_hwloc_base_cpu_set`) |
| `start`, `end`, `wallSeconds` | `Allrun` | Wall total de decomposePar + solver + reconstructPar |
| `status` | `Allrun` | `ok` si el log del solver termina en `End`, si no `failed` |
| `steps`, `simTime` | `Allrun` | Pasos de esta ejecución y tiempo físico inicial y final |
| `maxRSSkB` | `Allrun` | RAM pico del proceso MPI más grande, de `/usr/bin/time -v` |
| `coreUsPerCellStep` | `Allrun` | wall · np · 10⁶ / (nCells · steps) |
| `Cd` | `Allrun` | Último valor de `coefficient.dat` |

Las core-horas son `wallSeconds · np / 3600`. Para comparar con Nektar++:

- **Costo por grado de libertad:** `coreUsPerCellStep` contra el costo por
  DOF y por paso de Nektar++ (elementos · (P+1)³).
- **Costo por resultado:** core-horas hasta que Cd queda dentro de la
  tolerancia. Es la comparación que importa: un paso barato no garantiza
  un resultado barato.

Dos corridas que comparten la máquina se pisan el ancho de banda de memoria.
Para una comparación de tiempos, corré de a una o anotá qué más corría.

### Una corrida por socket

`newton` tiene 2 sockets de 14 núcleos físicos. Los CPUs lógicos 0–13 son
el socket 0 y 14–27 el socket 1; 28–55 son los hyperthreads. `runParallel`
no acepta opciones de `mpirun`, pero Open MPI lee las variables `OMPI_MCA_*`,
que pasan por `run.sh` y `Allrun`:

```bash
OMPI_MCA_hwloc_base_cpu_set=0-13  OMPI_MCA_hwloc_base_binding_policy=core ./run.sh sweep_m03.txt
OMPI_MCA_hwloc_base_cpu_set=14-27 OMPI_MCA_hwloc_base_binding_policy=core ./run.sh sweep_m06.txt
```

Cada corrida queda en su socket, con su memoria y su L3, un proceso por
núcleo físico. El ancho de banda deja de ser compartido entre corridas; la
frecuencia del turbo sí sigue dependiendo de la carga total.

`log.time.<fecha>` guarda la salida completa de `/usr/bin/time -v`.

## 3. Disco

Medido en `m03-coarse` (0.73 M celdas, 41.5 k pasos), antes de estos cambios:

| Qué | Antes | Ahora |
|---|---|---|
| `log.rhoCentralFoam` | 2.2 KB por paso | Sin la salida de `forceCoeffs`, y comprimido con `xz` al terminar (`log.<app>.<fecha>.xz`) |
| `residuals` | 290 B por paso | Una línea cada 50 pasos |
| `constant/polyMesh` | Copia por caso (unos 185 B por celda) | Enlace simbólico a `<meshes-dir>/<malla>`; copia sólo con `--refine` |
| Tiempos escritos | Los dos últimos, más las copias en `processor*` | El último (`purgeWrite 1`); `processor*` se borra si `reconstructPar` terminó |

Los logs comprimidos se leen con `xzgrep`, `xzless` o `xz -dc`.

> **CUIDADO** Con el enlace, borrar o rehacer una malla rompe todos los
> casos que la usan. `Allrefine` reemplaza el enlace por una copia antes de
> refinar, así no modifica la malla compartida.

Mientras corre, el log está sin comprimir y crece a unos 1–2 KB por paso.

## 4. Aviso al celular

`Allrun` manda un POST al terminar, con éxito o con error. Llega con el
nombre del caso, `ok` o `failed`, el tiempo físico, el Cd y las horas de
wall.

Con [ntfy](https://ntfy.sh): instalá la app, suscribite a un tema con un
nombre difícil de adivinar (el tema es público para quien lo conozca) y
escribí la URL en un archivo:

```bash
echo https://ntfy.sh/<tema> > ~/.config/aconcagua/notify
```

`Allrun` lee el archivo al terminar cada caso, así que sirve también para
corridas que ya están en marcha. La variable de entorno `NOTIFY_URL` tiene
prioridad sobre el archivo.

No avisa si se mata `run.sh` o el lanzador, o si se cae la máquina: en esos
casos `Allrun` no llega a su última parte.
