# Flujo de trabajo

De la geometría al Cd, en tres comandos.

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc

cd mesh && ./Allmesh && cd ..        # malla, ~5 min
./newCase.sh ~/runs/m08 --regime trans --Minf 0.8
cd ~/runs/m08 && ./Allrun 8
```

O el barrido entero de una: editá `sweep.txt` y corré `./run.sh`.

Poné los directorios de corrida **fuera de OneDrive y fuera de `/mnt/c`**: las
dos cosas hacen la E/S de OpenFOAM varias veces más lenta.

## 1. La malla

`mesh/Allmesh` corre cfMesh sobre `mesh/Aconcagua.stl` y deja
`mesh/constant/polyMesh`. Los seis pasos están en el script; los tamaños, en
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
aleta tiene 12 mm de espesor: con 10 capas la extrusión se enreda y aparecen
volúmenes negativos. Está medido en `docs/TROUBLESHOOTING.md`.

La malla actual son 862 k celdas, 95 % hexaedros, sin volúmenes negativos,
no-ortogonalidad máxima 78.6 y media 5.7.

Para mallar otra geometría: `./Allmesh mi_cohete.stl`. El STL tiene que traer
los solids nombrados `nosecone`, `body`, `boattail` y `fins` — `Allmesh` los
renombra a los patches que esperan los casos.

## 2. El caso

`newCase.sh` copia una plantilla más la malla, y escribe las condiciones de
vuelo que le pases.

```bash
./newCase.sh ~/runs/m02 --regime sub   --Uinf 68
./newCase.sh ~/runs/m12 --regime super --Minf 1.2 --np 16
./newCase.sh ~/runs/a05 --regime sub   --Uinf 68 --alpha 5 --refine tip
```

| opción | qué hace |
|---|---|
| `--regime sub\|trans\|super` | qué plantilla. Ver [SOLVERS.md](SOLVERS.md) |
| `--np N` | descomposición |
| `--refine tip\|fins` | refinamiento local en las aletas, ver sección 3 |
| condiciones | `--alpha --beta --Ti --nuRatio`; `sub`: `--Uinf --nu --rhoInf`; `trans`/`super`: `--Minf --pInf --Tinf` |

El caso nunca se corre en la plantilla: `newCase.sh` se niega a pisar un
directorio que ya existe.

## 3. El barrido

`sweep.txt` es una línea por corrida:

```
m02  --regime sub    --Uinf 68
m08  --regime trans  --Minf 0.8
m18  --regime super  --Minf 1.8
```

`./run.sh` las arma y corre en serie bajo `runs/`. Saltea las que ya existen,
así que si se corta, volvés a lanzarlo y sigue donde estaba.

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

Secuencia: `restore0Dir` → `decomposePar -force` → `potentialFoam` (campo
inicial, ayuda mucho a la convergencia) → `simpleFoam` → `reconstructPar
-latestTime`. Logs en `log.<aplicación>`.

`controlDict`: `endTime 3000` iteraciones, `residualControl` en `fvSolution`
corta antes (p 1e-5, U/k/omega 1e-6). `startFrom latestTime`, así `./Allrun`
sobre una corrida existente continúa en vez de empezar de cero (borrar `[1-9]*`
y `processor*` para reiniciar, o `./Allclean` para volver a la plantilla).

Orden de magnitud: la malla fina de un cuarto (3.5 M) a 8 procesos hace
~1 iteración/s; 2000 iteraciones son ~40 min. La completa fina (14 M) es
4× eso y necesita ~30 GB de RAM entre solver y reconstrucción.

## 6. Resultados

| Qué | Dónde |
|---|---|
| `Cd`, `Cl`, `Cm` por iteración | `postProcessing/forceCoeffs1/0/coefficient.dat` |
| fuerzas en N (presión, viscosa) | `postProcessing/forces1/0/force.dat`, `moment.dat` |
| y+ por patch | `log.simpleFoam` (`grep -A5 "yPlus yPlus write"`) y campo `yPlus` en cada tiempo escrito |
| `Cp`, `wallShearStress` | campos en los tiempos escritos |
| residuales | `postProcessing/residuals/0/solverInfo.dat` |
| ParaView | `touch caso.foam; paraview caso.foam` (o `paraFoam`) |

## 7. Estudios típicos

**Barrido de Mach.** Es para lo que existe `sweep.txt`. Una malla por punto:

```bash
./run.sh
```

**Convergencia de malla.** Subí o bajá `maxCellSize` en `mesh/system/meshDict`
por factores de 2, dejando los niveles quietos: eso escala todo el campo sin
tocar la resolución relativa. Corré el mismo punto de Mach en cada malla y
graficá Cd contra el tamaño de celda.

**Barrido de ángulo de ataque.**

```bash
for a in 0 2 4 6 8; do
    ./newCase.sh ~/runs/a$a --regime sub --Uinf 68 --alpha $a --np 8
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

- `mesh/constant/polyMesh` no se versiona: se reconstruye con `./Allmesh` en
  cinco minutos. Lo que se versiona es el STL y el `meshDict`.
- Los directorios de corrida van fuera de la repo, fuera de OneDrive y fuera
  de `/mnt/c`.
- Un cambio en el `meshDict` que mueva la malla se commitea junto con el
  número de celdas y la salida de `checkMesh` que produjo.

## Chuleta

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc

cd mesh && ./Allmesh && cd ..              # malla
./newCase.sh ~/runs/NOMBRE --regime sub --Uinf 68 --np 8
cd ~/runs/NOMBRE && ./Allrun 8

./run.sh                                   # el barrido entero de sweep.txt

./Allrefine tip                            # refinamiento local, antes de Allrun
./Allclean                                 # volver a la plantilla
```
