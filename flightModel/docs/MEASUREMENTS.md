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
| `meshHash`, `flow` | `newCase.sh` | Huella de la malla y definición física del caso (§6) |
| `host`, `cpu`, `np`, `binding` | `Allrun` | Máquina, procesos MPI y socket fijado. Con un solo socket se fija el 0 aunque no haya `--socket`, si `np` entra en los núcleos físicos; `none` en el archivo `socket` lo desactiva (para una prueba al lado de otra corrida) |
| `start`, `end`, `wallSeconds` | `Allrun` | Wall total de decomposePar + solver + reconstructPar |
| `solverWallSeconds` | `Allrun` | Wall del solver solo |
| `status` | `Allrun` | `ok` si el log del solver termina en `End`, si no `failed` |
| `steps`, `simTime` | `Allrun` | Pasos de esta ejecución y tiempo físico inicial y final |
| `maxRSSkB` | `Allrun` | RAM pico del proceso MPI más grande, de `/usr/bin/time -v` |
| `coreUsPerCellStep` | `Allrun` | wall · np · 10⁶ / (nCells · steps) |
| `cyclesPerCellStep` | `Allrun` | Ciclos de CPU de todos los procesos / (nCells · steps) (§5) |
| `Cd` | `Allrun` | Último valor de `coefficient.dat` |
| `machine`, `measured` | `Allrun` | Estado de la máquina al arrancar y lo medido durante la corrida (§5) |

Las core-horas son `wallSeconds · np / 3600`. Para comparar con Nektar++
están las §5 y §6.

Dos corridas que comparten la máquina se pisan el ancho de banda de memoria.
Para una comparación de tiempos, corré de a una o anotá qué más corría.

### Una corrida por socket

`newton` tiene 2 sockets de 14 núcleos físicos; los CPUs lógicos 28–55 son
hyperthreads. `--socket <S>` en la línea del sweep (o en `newCase.sh`) fija
la corrida al socket `S`:

```
m03-h1  --np 8 --socket 0 --mesh h1 ...
m06-h1  --np 8 --socket 1 --mesh h1 ...
```

`newCase.sh` escribe `S` en el archivo `socket` del caso. `Allrun` arma un
`rankfile` con el proceso `i` en el núcleo físico `i` de ese socket, y se lo
pasa a Open MPI por `OMPI_MCA_rmaps_rank_file_path`. Llama a `mpirun`
directamente, con `--bind-to core`, y no a `runParallel`: `runParallel`
siempre agrega `--oversubscribe`, y con eso Open MPI 5 ignora el `rankfile` y
deja cada proceso suelto en todos los CPUs (nemepad, 2026-10-17). Se detiene
si `np` supera los núcleos del socket. Para comprobar el binding:
`taskset -cp <pid>` de cada proceso tiene que dar los dos hilos de un solo
núcleo.

Cada corrida queda con su memoria, su L3 y un núcleo físico por proceso, sin
hyperthreads ni migraciones. La frecuencia del turbo sí sigue dependiendo de
la carga total.

`--cpu-set` (`OMPI_MCA_hwloc_base_cpu_set`) no alcanza: Open MPI 4.1 deja
cada proceso suelto en todo el socket aunque se pida `--bind-to core`.

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

No avisa si se mata `runCase.sh`, o si se cae la máquina: en esos
casos `Allrun` no llega a su última parte.

## 5. Medición de la máquina

El wall de una corrida depende de la máquina tanto como del solver: perfil
de energía, turbo, temperatura, suspensión. Para comparar OpenFOAM con
Nektar++ hace falta saber a qué velocidad corrió cada uno. `Allrun` lo
mide sin intervenir en el solver.

### Qué corre al lado del solver

| Archivo | Herramienta | Contenido |
|---|---|---|
| `log.perf.<fecha>.xz` | `perf stat -I 30000` sobre el solver y sus procesos MPI | Cada 30 s: `task-clock`, `cycles`, `instructions`, `stalled-cycles-frontend`, `cache-misses`, migraciones y cambios de contexto |
| `log.machine.<fecha>` | `machineLog.py sample`, cada 30 s | Reloj de pared, `CLOCK_MONOTONIC`, `CLOCK_BOOTTIME`, MHz medio, mínimo y máximo de los núcleos (`cpuinfo_avg_freq`), carga, temperatura de CPU, de cada DIMM y del NVMe, ventiladores y `platform_profile` |
| `log.energy.<fecha>.xz` | `perf stat -a -e power/energy-pkg/` | Energía del paquete en J, cada 30 s. Sólo si `perf` puede contar en todo el sistema |

### Qué queda en `run.info`

`machine` es el estado al arrancar: kernel, sockets, núcleos, hilos por
núcleo, RAM, governor, EPP, boost, `scaling_max_freq`, `platform_profile`,
nivel del ventilador, `perf_event_paranoid`, versión y compilación de
OpenFOAM, versión de Open MPI.

`measured` es lo que pasó durante el solver:

| Entrada | Qué |
|---|---|
| `suspendedSeconds` | Δ`BOOTTIME` − Δ`MONOTONIC`: el tiempo que la máquina estuvo suspendida. Si no es 0, el wall no sirve como costo |
| `mhzMean`, `mhzMinOfMean` | Reloj medio de todos los CPUs lógicos, y el peor promedio de 30 s |
| `tCpuMax`, `tDimm<i>Max`, `tNvmeMax` | Temperaturas máximas. La DDR5 dobla el refresh arriba de 85 °C |
| `profiles` | Los `platform_profile` que aparecieron. Si no es sólo `performance`, la corrida no compara |
| `cpuSeconds` | CPU de todos los procesos del solver |
| `cycles`, `instructions` | Totales de todos los procesos |
| `GHzEffective` | `cycles / cpuSeconds`: el reloj real al que corrió el solver, con turbo y throttling incluidos |
| `IPC` | `instructions / cycles`. Bajo (≲ 1) en un código limitado por memoria |
| `frontendStallFraction` | Ciclos con el frontend parado / ciclos |
| `cacheMisses`, `cpuMigrations`, `contextSwitches` | Para detectar ruido: migraciones altas = binding roto |
| `energyPkgJ` | Energía del paquete durante el solver |

`cyclesPerCellStep` (§2) es el costo que no depende del reloj, el que se
pone al lado del costo por DOF de Nektar++. Incluye el arranque (leer la
malla), que en una corrida de miles de pasos es despreciable.

### Requisitos

- `kernel.perf_event_paranoid ≤ 1` para contar ciclos de los propios
  procesos. Sin eso `Allrun` avisa con `!!` y corre igual, sin `cycles`.
- La energía necesita `≤ 0` (o `CAP_PERFMON`): `perf` la lee de RAPL para
  todo el paquete. Con `1` no queda `energyPkgJ`.

  ```bash
  sudo sysctl kernel.perf_event_paranoid=0
  echo kernel.perf_event_paranoid=0 | sudo tee /etc/sysctl.d/90-perf.conf
  ```

### Protocolo para que el costo compare

1. `platform_profile` en `performance`. `Allrun` avisa si no lo está.
2. Ventiladores fijos (en una ThinkPad, `level 7` con
   `thinkpad_acpi fan_control=1`). Recargar el módulo `thinkpad_acpi`
   puede dejar el perfil en `custom`: volvé a ponerlo en `performance`.
3. Que la tapa no suspenda (`HandleLidSwitch=ignore` en logind, o un
   `systemd-inhibit --what=handle-lid-switch:sleep`).
4. Una corrida por vez y nada pesado al lado: dos corridas se pisan el
   ancho de banda de memoria.
5. El mismo `np` y el mismo binding en los dos solvers: un proceso por núcleo
   físico, sin hyperthreads.

Lo que el protocolo no controla queda registrado en `measured`. Antes de
usar el costo de una corrida, mirá `suspendedSeconds`, `profiles` y
`GHzEffective`.

### No dividir por los GHz

Corregir el wall por el reloj supone que el tiempo escala con 1/f. Vale para
un código limitado por cálculo, no para uno limitado por memoria, y la
diferencia entre los dos solvers es justamente esa: OpenFOAM, de bajo orden,
depende más del ancho de banda; Nektar++ a P alto, más del reloj. La forma de
medirlo es un barrido de reloj (`scaling_max_freq` en 2, 3 y 4 GHz, el mismo
caso corto en cada solver) y ajustar t(f) = t_mem + C/f. Eso no está
automatizado todavía.

### Corridas sin esta medición

Las anteriores a este cambio tienen sólo `timeInfo.dat`: CPU y wall del
proceso 0 cada 50 pasos. Un intervalo con wall ≫ CPU es una suspensión; uno
con CPU muy por encima de la mediana, reloj bajo o carga al lado. Sacando
esos intervalos, el costo por paso de las corridas h0 en nemepad da
20–22 µs·núcleo por celda y paso (2026-10-17). No hay ciclos ni GHz.

## 6. Comparar con Nektar++

"El mismo caso" es la misma malla, las mismas condiciones de entorno y el
mismo campo inicial. Los esquemas pueden diferir; lo que sigue es lo que
cada corrida de OpenFOAM deja para que la de Nektar++ lo reproduzca, y lo
que hay que medir igual en las dos.

### Lo que tiene que ser igual

| Qué | Dónde está | Nota |
|---|---|---|
| Malla | `meshHash` en `run.info`; `<malla>/meshHash` | sha256 (16 caracteres) de `points`, `faces`, `owner`, `neighbour`, `boundary`. La malla de Nektar++ sale de la misma `polyMesh` |
| Media malla | `meshInfo`: `sector half` | Plano de simetría z = 0. `Aref` es la mitad del área frontal |
| Corriente libre | `flow` en `run.info` | `Minf`, `alpha`, `pInf`, `Tinf`, `Uinf`, `rhoInf`, `muInf`, las componentes de U |
| Reynolds | `RePerMetre`, `ReD` (diámetro), `ReL` (largo del cuerpo) | Con `rhoInf`, `Uinf`, `muInf` |
| Gas | `Rgas`, `gammaGas`, `Cp`, `RgasThermo` | `flowDerived` usa R = 287.05 y el termo, 8314.46/28.96 = 287.10: 0.02 % de diferencia en `rhoInf` |
| Viscosidad | `sutherland (As Ts)` | μ = As·√T/(1 + Ts/T) |
| Conductividad | `PrEucken` | OpenFOAM con `sutherland` saca κ con la corrección de Eucken: Pr ≈ 0.69, no 0.72. Nektar++ tiene que usar ese Pr |
| Turbulencia | `turbulence`, `Prt`, `Ti`, `nuRatio`, `kInlet`, `omegaInlet` | k-ω SST con `nutUSpaldingWallFunction` en las paredes |
| Pared | `wallT`, `wallNut` | Adiabática (`zeroGradient` en T), no deslizamiento |
| Campo inicial | `initial` | `uniform free stream` o la corrida de la que se mapeó |
| Referencias | `FOcommon`: `CofR (0 0 0)` (punta de la ojiva), `dragDir` y `liftDir` en ejes viento, `lRef`, `Aref` | `Cd` y `Cl` en ejes viento |

### Lo que se compara

| Magnitud | Archivo en OpenFOAM |
|---|---|
| Historia de `Cd`, `Cl`, `CmPitch`, cada 50 pasos | `forceCoeffs1/*/coefficient.dat` |
| Fuerza de presión y viscosa | `forces1/*/force.dat` |
| Cp, τw e y+ a lo largo del cuerpo | `samples/*/skinZ0_*.csv` (corte 0.5 mm sobre z = 0, con x, y, z: barlovento y < 0, sotavento y > 0) y `skinY0_*.csv` (sobre y = 0, de costado). Los dos planos pasan entre las aletas |
| Perfiles de estela | `samples/*/wake05`, `wake1`, `wake2`, `wake5`: líneas en y de −0.6 a 0.6 m sobre z = 0, a 0.5, 1, 2 y 5 diámetros detrás de la base |
| Eje | `samples/*/axis_*.csv`: de x = −0.5 m a 20 diámetros detrás de la base |
| Control de corriente libre | `samples/*/freeStream_*.csv`: un punto 1 m delante de la punta |
| Piel y plano de simetría completos | `surfaces/*/` (VTK) |
| Costo | `run.info`: `cyclesPerCellStep`, `cpuSeconds`, `energyPkgJ`, `GHzEffective` |

Las líneas se escriben al final (`onEnd`), con `cellPoint`, en CSV.

### Cómo medir igual en los dos

- **Grados de libertad.** OpenFOAM: nCells por variable. Nektar++:
  elementos · (P+1)³ por variable en hexaedros (otro conteo en prismas y
  tetraedros). El costo se compara por DOF y por paso, en ciclos.
- **Convergencia.** Cada solver corta con su propio criterio y su propio
  paso. Para el costo hasta el resultado, aplicá el mismo criterio a las dos
  historias de `Cd` en el post, contra el costo acumulado, y no uses los
  residuos: se normalizan distinto en cada código.
- **Refinamiento.** Nektar++ suma refinamiento p a las mallas h1–h3. Cada
  par (h, P) es un caso; el costo por DOF permite ponerlos en la misma
  curva que h1–h3 de OpenFOAM.

### Arranque desde otra corrida

`Allrun` corta si `mapFields` o `changeDictionary` no terminan en `End`.
Antes sólo buscaba `FOAM FATAL`, y un segfault de `mapFields` dejaba
arrancar desde el flujo uniforme diciendo lo contrario. Pasa al mapear
sobre la misma malla h0 → h0 (al interpolar `Ma`, 2026-10-17), que es lo que
hace la cadena supersónica pimple → central.
