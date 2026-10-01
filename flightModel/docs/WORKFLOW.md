# Flujo de trabajo

De la geometría al Cd, en tres comandos.

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc

mesh/Allmesh ~/meshes/base          # malla, ~2 min
./newCase.sh ~/runs/m08 --mesh ~/meshes/base --Minf 0.8
cd ~/runs/m08 && ./Allrun 8
```

O el barrido entero de una: editá `sweep.txt` y corré `./run.sh`.

Poné los directorios de corrida **fuera de OneDrive y fuera de `/mnt/c`**: las
dos cosas hacen la E/S de OpenFOAM varias veces más lenta.

## 1. La malla

`mesh/Allmesh <dir>` corre cfMesh sobre `mesh/stl/Aconcagua.stl` y deja la
malla en `<dir>/constant/polyMesh`, junto con una copia de `system/` y
`meshInfo`: cada malla guarda el `meshDict` que la hizo, así que pueden
convivir varias (base, fina, otra geometría). Nunca pisa un directorio que
existe. Los seis pasos están en el script; los tamaños, en
`mesh/system/meshDict`.

Todo en `meshDict` es un **nivel**: un entero que parte el tamaño de celda al
medio, a partir de un único `maxCellSize`.

| nivel | tamaño |
|---|---|
| 4 | 100 mm |
| 6 | 25 mm |
| 8 | 6.25 mm |
| 9 | 3.13 mm |
| 10 | 1.56 mm |

Tres grupos de perillas y nada más:

- `localRefinement` — un nivel y un espesor de banda por patch
  (`cone`, `walls`, `tail`, `fins`).
- `objectRefinements` — cajas de refinamiento de volumen.
- `boundaryLayers` — `nLayers` y `thicknessRatio`, y `nLayers` por patch.

**Las aletas van con `nLayers 3`, no más.** El borde de ataque es un filo y la
aleta tiene 6 mm de espesor: con 10 capas la extrusión se enreda y aparecen
volúmenes negativos. Está medido en `docs/TROUBLESHOOTING.md`.

La malla actual son 2.1 M celdas, 96 % hexaedros, sin volúmenes negativos.
Detrás de la base la estela baja de a un nivel en tres cilindros coaxiales
(`baseWake` 6.25 mm, `nearWake` 12.5 mm, `midWake` 25 mm) antes de entrar a
`wake`.

Para mallar otra geometría: ponela en `mesh/stl/` y `mesh/Allmesh ~/meshes/otro mesh/stl/mi_cohete.stl`. El STL tiene que traer
los solids nombrados `nosecone`, `body`, `boattail` y `fins` — `Allmesh` los
renombra a los patches que esperan los casos.

## 2. El caso

`newCase.sh` copia una plantilla más la malla, y escribe las condiciones de
vuelo que le pases.

```bash
./newCase.sh ~/runs/m03 --mesh ~/meshes/base --Minf 0.3
./newCase.sh ~/runs/m12 --mesh ~/meshes/base --Minf 1.2 --np 16
./newCase.sh ~/runs/a05 --mesh ~/meshes/base --Minf 0.8 --alpha 5 --refine tip
```

| opción | qué hace |
|---|---|
| `--np N` | descomposición |
| `--refine tip\|fins` | refinamiento local en las aletas, ver sección 3 |
| malla | `--mesh <dir>`, obligatorio |
| condiciones | `--Minf --pInf --Tinf --alpha --beta --Ti --nuRatio` |

El caso nunca se corre en la plantilla: `newCase.sh` se niega a pisar un
directorio que ya existe.

Si la plantilla trae archivos `*.j2`, `newCase.sh` los renderiza con
`jinja2` contra su `config.json`. Es la rampa de Courant: [SOLVERS.md
§2.3](SOLVERS.md).

## 3. El barrido

`sweep.txt` es una línea por corrida:

```
m03     --mesh base --Minf 0.3
m08     --mesh base --Minf 0.8
m03-h2  --mesh h2 --maxCellSize 0.8 --Minf 0.3
```

`./run.sh [sweep.txt] [runs] [mallas]` las arma y corre en serie bajo `runs/`.
Saltea las que ya existen, así que si se corta, volvés a lanzarlo y sigue
donde estaba.

`--mesh` es obligatorio. Un nombre suelto vive en `~/meshes/`. Si la malla no
existe, `run.sh` la construye antes con `mesh/Allmesh`, aplicando
`--maxCellSize` y `--finLevel` (nivel de `fins`) a la copia del `meshDict`; las
líneas siguientes que la nombran la reusan. Si una línea pide una malla que ya
existe con otro `maxCellSize` o `finLevel`, el barrido se detiene en vez de
correr sobre la malla equivocada. `Allmesh` también se detiene si la malla sale
con volúmenes negativos.

**Una malla por punto de Mach.** `y1` sale de la extrusión de capas contra la
celda de superficie, no de un y+ objetivo, así que el y+ real cambia con la
velocidad. Medilo en cada corrida (`docs/TROUBLESHOOTING.md`, sección y+) y
ajustá `nLayers` si hace falta.

## 4. Condiciones de vuelo

Un archivo: `system/flowConditions`. Números planos.

```
Uinf   100;    // m/s
alpha  0;      // deg, ángulo de ataque: U rota de +x hacia +y. Requiere half o full
beta   0;      // deg, deslizamiento: hacia +z. Requiere full
nu     1.5e-05;
rhoInf 1.225;
Ti     0.005;  // intensidad de turbulencia
nuRatio 5;     // nut/nu en la corriente libre
```

`system/flowDerived` convierte eso en `(Ux Uy Uz)`, direcciones de drag y
lift, `k`, `omega`, `nut` con `#eval`; los campos de `0.orig/`,
`transportProperties` y `controlDict` incluyen `flowDerived`. No hace falta
tocar ningún campo para cambiar la velocidad o el ángulo.

Editarlo con `foamDictionary system/flowConditions -entry alpha -set 5` o con
un editor. **No** usar `foamDictionary -set` sobre `controlDict` ni sobre los
campos: reescribe el archivo con los `#include` expandidos y los `#eval`
evaluados, y las condiciones dejan de seguir a `flowConditions`.

`Allrun` comprueba `alpha`/`beta` contra el `sector` del `meshInfo` y se
niega a arrancar si la malla no puede representar ese ángulo.

## 5. Correr

```bash
./Allrun          # paralelo, numberOfSubdomains de decomposeParDict (8)
./Allrun 16       # 16 procesos (actualiza decomposeParDict)
./Allrun 1        # serie
```

Secuencia: `restore0Dir` → `decomposePar -force` → `rhoCentralFoam` →
`reconstructPar -latestTime`. Logs en `log.<aplicación>`.

`endTime` es tiempo físico; el corte real es el `runTimeControl` sobre el
promedio de `Cd`, que se activa en `minIter` ([SOLVERS.md §2.3](SOLVERS.md)).
`startFrom latestTime`, así `./Allrun` sobre una corrida existente continúa en
vez de empezar de cero (borrar `[1-9]*` y `processor*` para reiniciar, o
`./Allclean` para volver a la plantilla).

Orden de magnitud: la malla completa de 2.3 M a 8 procesos hace ~2.3 s por
paso.

## 6. Resultados

| Qué | Dónde |
|---|---|
| `Cd`, `Cl`, `Cm` por iteración | `postProcessing/forceCoeffs1/0/coefficient.dat` |
| fuerzas en N (presión, viscosa) | `postProcessing/forces1/0/force.dat`, `moment.dat` |
| y+ por patch | `log.rhoCentralFoam` (`grep -A5 "yPlus yPlus write"`) y campo `yPlus` en cada tiempo escrito |
| `Cp`, `wallShearStress` | campos en los tiempos escritos |
| residuales | `postProcessing/residuals/0/solverInfo.dat` |
| ParaView | `touch caso.foam; paraview caso.foam` (o `paraFoam`) |

## 7. Estudios típicos

**Barrido de Mach.** Es para lo que existe `sweep.txt`. Una malla por punto:

```bash
./run.sh
```

**Convergencia de malla.** Se varía sólo `maxCellSize`, con una razón
constante, y los niveles quietos: eso escala todo el campo sin tocar la resolución
relativa, que es lo que pide una extrapolación de Richardson. En `sweep.txt`:

```
m03-h1  --mesh h1 --maxCellSize 1.6 --Minf 0.3
m03-h2  --mesh h2 --maxCellSize 1.13 --Minf 0.3
m03-h3  --mesh h3 --maxCellSize 0.8 --Minf 0.3
```

La razón es √2 (el GCI de Celik et al. 2008 pide más de 1.3): cada paso
multiplica las celdas por hasta 2√2 ≈ 2.8, así que `h3` queda en ~17 M como
mucho. Con razón 2 serían ×8 por paso, y la tercera malla pasaría los 100 M.
El primer espesor de capa se achica con la celda de superficie, así que el y+
también cambia entre mallas: medilo en las tres. El OpenFOAM ToolChain
hace sus mallas R1–R6 a mano, cambiando a la vez `maxCellSize` y la cantidad
de cajas; sirven para elegir una malla, pero la razón de refinamiento no es
constante y no alcanzan para estimar el error de discretización.

**Barrido de ángulo de ataque.**

```bash
for a in 0 2 4 6 8; do
    ./newCase.sh ~/runs/a$a --mesh ~/meshes/base --Minf 0.8 --alpha $a --np 8
    ( cd ~/runs/a$a && ./Allrun 8 )
done
```

**Refinar más las aletas.** Subí `fins` a nivel 10 en `localRefinement`, o usá
`./Allrefine tip` en el caso ya armado, que refina 2:1 una banda alrededor de
la punta sin volver a mallar.

**Estela transitoria.** Cambiá el solver a `pimpleFoam` (`ddtSchemes backward`,
`deltaT` para Co ≈ 1 en la estela cercana: con celdas de 6 mm y 100 m/s,
~5e-5 s), y agrandá la caja `wake` del `meshDict`.

## 8. Convenciones

- Las mallas no se versionan y viven fuera de la repo: se reconstruyen con
  `mesh/Allmesh <dir>` en
  cinco minutos. Lo que se versiona es el STL y el `meshDict`.
- Los directorios de corrida van fuera de la repo, fuera de OneDrive y fuera
  de `/mnt/c`.
- Un cambio en el `meshDict` que mueva la malla se commitea junto con el
  número de celdas y la salida de `checkMesh` que produjo.

## Chuleta

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc

mesh/Allmesh ~/meshes/base                 # malla
./newCase.sh ~/runs/NOMBRE --mesh ~/meshes/base --Minf 0.8 --np 8
cd ~/runs/NOMBRE && ./Allrun 8

./run.sh                                   # el barrido entero de sweep.txt

./Allrefine tip                            # refinamiento local, antes de Allrun
./Allclean                                 # volver a la plantilla
```
