# Flujo de trabajo

De la geometría al `Cd`, en orden. Cada paso es un comando; cada archivo que
se edita está nombrado. Si venís de OpenFOAM: la malla reemplaza a
`blockMesh` + `snappyHexMesh`, el caso es un `simpleFoam` estándar con
`kOmegaSST` y wall functions, y los `Allmesh`/`Allrun`/`Allclean` hacen lo
que hacen en cualquier tutorial.

```
aconcaguaGeom.py     geometría de registro (5 números + aletas)
       │
meshParams.py  ─┐
presets/*.py    ├─►  python build.py [--preset] [--scale] [--sector] [--set]
línea de comando┘           │
                            ├─► output/<nombre>.msh         malla gmsh v2.2 ASCII
                            ├─► output/<nombre>.meshInfo    dict OpenFOAM: sector, Aref, patches
                            └─► output/<nombre>.params.py   parámetros usados (reproducible)
                                       │
./newCase.sh <dir> <msh> [--regime R]  │   copia case-<R>/ y corre Allmesh
                                       ▼
                     <dir>/  Allmesh  →  gmshToFoam + fixPatchTypes + checkMesh + renumberMesh
                             Allrefine [tip|fins]    (opcional) topoSet + refineMesh
                             system/flowConditions   (Uinf, alpha, beta, nu, ...)
                             Allrun [N]  →  potentialFoam + simpleFoam (paralelo)
                             postProcessing/forceCoeffs1/0/coefficient.dat
```

## 0. Entorno

**Windows + WSL** (la configuración del equipo): la malla se construye con
el Python de Windows (`python`, con `gmsh` instalado por pip) o con el de WSL
(`python3`); OpenFOAM sólo existe en WSL. Los scripts de shell (`Allmesh`,
`Allrun`, `newCase.sh`) se corren desde WSL con OpenFOAM cargado:

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
```

El repo se ve desde WSL en `/mnt/c/Users/<usuario>/OneDrive/CFD-ITBA-ROCKETRY`.
Construir la malla ahí está bien; **correr OpenFOAM ahí no**: `/mnt/c` es
varias veces más lento para I/O y OneDrive intenta sincronizar cada archivo
de `processor*/`. Las corridas van a un directorio de WSL (`~/runs/...`).

**Linux nativo**: todo en la misma shell.

## 1. Construir la malla

```bash
cd mesh
python build.py --plan                      # el plan: celdas, y1, y+, capas en la BL
python build.py --preset coarse             # construye, audita y escribe
```

`--plan` imprime el plan de celdas resuelto (radial, axial, azimutal y total
estimado) sin construir nada. Usalo siempre antes de una malla grande: tarda
un segundo y el total estimado difiere del real en menos de 1 %.

### Presets

Un preset es un archivo `presets/<nombre>.py` con líneas `PARAMETRO = valor`.
Se apilan (`--preset medium --preset fintip`) y se combinan con los flags.

| Preset | Celdas / cuadrante | Tiempo build | Para qué |
|---|---|---|---|
| `smoke` | 31 k | 25 s | probar el pipeline; todo escalado ×6, y+ ≈ 200, aletas de una celda: no sirve para resultados |
| `coarse` | 530 k | 45 s | primera mirada al campo. **Mantiene y+ = 32 y la aleta a resolución completa**; el resto ×3 |
| `medium` | 1.22 M | 1.5 min | resolución de trabajo; y+ y aleta como `fine`, el resto ×2 |
| `fine` | 3.54 M | 2 min | malla de referencia (los defaults de `meshParams.py`) |
| `fintip` | +36 % | | 6 mm radiales en la punta de las aletas (`FIN_H_R`) |
| `wake_unsteady` | +65 % | | estela cercana más larga y fina, para `pimpleFoam` |

`coarse` y `medium` llevan `H_SCALE_WALL = False` y `H_SCALE_FIN = False`:
la primera celda de cada pared, y la aleta en sus tres direcciones (5 mm sobre
la cuerda, 0.5 mm normal, 12 mm radiales en la punta), no se engrosan. Una
aleta de 12 mm de espesor engrosada ×3 queda de una celda: se ve como un bulto
en ParaView y su patch `fins` es 30 % más grande que la planform.

Tiempos en un escritorio de 16 núcleos; el build es serie (gmsh) y escala
con el número de celdas. Media malla ×2, completa ×4.

### Flags

| Flag | Efecto |
|---|---|
| `--scale H` | `H_SCALE`: multiplica todo tamaño de celda por H (celdas ≈ 1/H³). `--scale 2` es 1/8 de las celdas. |
| `--sector quarter\|half\|full` | 90°, 180° o 360°. Ver [SECTORS_AND_AOA.md](SECTORS_AND_AOA.md). |
| `--fins` / `--no-fins` | con o sin aletas (cuerpo de revolución limpio) |
| `--set KEY=VALUE` | cualquier parámetro; `VALUE` es un literal de Python: `--set ZONE_R=[0.28,0.35,2,11.775]` |
| `--out FILE.msh` / `--name NOMBRE` | dónde y cómo se llama la salida. Default `output/aconcagua_<sector>_s<scale>[_nofins].msh` |
| `--no-write` | construye y audita, no escribe (para probar parámetros) |
| `--no-audit` | omite la auditoría de la malla ensamblada (mitad/completa); el cuadrante se audita siempre |

Precedencia, de menor a mayor: defaults de `meshParams.py` < `--preset` <
`--scale`/`--sector`/`--fins` < `--set`. Un nombre de parámetro mal escrito es
un error, no un no-op silencioso.

### Qué mirar en la salida

```
  wall projection (max |r - r_exact|)        <- después del snap debe ser 0.00 nm
  fins: 95,616 nodes deformed ...            <- las aletas se introducen deformando theta
  ---- quadrant audit ----
  faces in patches      266,244   OK        <- toda cara de borde tiene patch (nada en defaultFaces)
  a face shared by 3+ cells   False
  cells with volume <= 0          0
  non-orthogonality  max   75.85 deg         <- medido como lo mide checkMesh, ver nota
  skewness  max  2.75  (internal 0.7, boundary 2.75; > 4: 0)   <- las mismas definiciones que checkMesh
  fin: 5,452 faces, wetted area per half-fin 401.80 cm2 vs 376.01 exact (+6.9 %)
  wrote .../aconcagua_quarter_s1.msh  (523 MB)
```

`build.py` se niega a escribir si la auditoría encuentra celdas de volumen
negativo, caras compartidas por tres celdas o caras de borde sin patch. Para
mallas de mitad o completas, la auditoría del ensamblado es la prueba de que
la costura entre cuadrantes cerró: si un par de caras no se fusionó,
aparecería como cara de borde sin patch.

El área mojada de la aleta sale un poco más grande que la exacta porque el
contorno del patch está cuantizado a una celda. Es un sesgo conocido, medido y
convergente; qué significa para los coeficientes está en la sección 6.

**Nota sobre la auditoría.** El skewness usa las definiciones de OpenFOAM
(caras internas y de borde) y coincide con `checkMesh` a la segunda cifra. La
no-ortogonalidad no: la auditoría usa centros de celda promedio de vértices y
`checkMesh` centroides ponderados por volumen. En la malla fina: 75.9° máx.
(auditoría) contra 86.7° máx. (`checkMesh`), con 7.5 k caras > 70° de 10.7 M.
Los números de `checkMesh` son los que valen y están en
[MESH_DESIGN.md](MESH_DESIGN.md#current-quality).

### Los tres archivos de salida

- `<nombre>.msh` — 520 MB para el cuarto fino. No se versiona, no se copia:
  se regenera.
- `<nombre>.meshInfo` — diccionario OpenFOAM con `sector`, `Aref`, `lRef`,
  `wallPatches`, `patchTypes`, `nCells`, `y1`. `Allmesh` lo copia a
  `constant/meshInfo` y `controlDict` lo incluye: el `forceCoeffs` toma de
  ahí el `Aref` correcto para el sector y la lista de patches de pared.
- `<nombre>.params.py` — todos los parámetros con los que se construyó, en
  formato preset. `python build.py --preset output/<nombre>.params.py`
  reconstruye exactamente esa malla. Guardalo junto con los resultados.

## 2. Del `.msh` al caso

Hay **tres plantillas**, una por régimen de velocidad, y ninguna se corre en
el repo. `newCase.sh` copia la que elijas, fija las condiciones de vuelo y
convierte la malla:

| `--regime` | Directorio | Solver | Rango |
|---|---|---|---|
| `sub` (default) | `case-subsonic/` | `simpleFoam` | M < 0.3 |
| `trans` | `case-transonic/` | `rhoSimpleFoam` con `transonic yes` | 0.3 a 1.2 |
| `super` | `case-supersonic/` | `rhoCentralFoam` | M > 1.2 |

Los modelos de cada una, con sus referencias, están en
[SOLVERS.md](SOLVERS.md). Las dos compresibles **no están validadas**: corren,
pero sus números no tienen respaldo hasta hacer [VALIDATION.md](VALIDATION.md).

En las compresibles la condición se da como número de Mach y la velocidad se
deriva, así que `--Uinf` no existe ahí y `--Minf` no existe en la subsónica.
`newCase.sh` rechaza la combinación equivocada en vez de ignorarla.

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
./newCase.sh ~/runs/a05  mesh/output/aconcagua_half_s2.msh  --alpha 5 --np 8
```

Opciones: `--alpha --beta --Uinf --nu --rhoInf --Ti --nuRatio` (van a
`system/flowConditions`) y `--np` (`decomposeParDict`). El script nunca
sobreescribe un directorio existente.

A mano es lo mismo:

```bash
cp -r case-subsonic ~/runs/a05 && cd ~/runs/a05
./Allmesh /ruta/a/aconcagua_half_s2.msh
foamDictionary system/flowConditions -entry alpha -set 5
```

`Allmesh` hace, en orden: `gmshToFoam` → copia el `.meshInfo` a
`constant/meshInfo` → `fixPatchTypes.py` → `checkMesh -allGeometry -allTopology`
→ `renumberMesh -overwrite`.

`fixPatchTypes.py` **no es opcional**: `gmshToFoam` crea todos los patches
como `type patch`. `symm` como `patch` falla ruidosamente; `cone`/`walls`/
`tail`/`fins` como `patch` en vez de `wall` **no falla**: las wall functions,
`yPlus` y `forceCoeffs` dan resultados silenciosamente incorrectos. Los tipos
salen del `patchTypes` de `constant/meshInfo`, así que una malla completa (sin
`symm`) se maneja sola.

### Qué es normal en `checkMesh`

Con `-allGeometry`, la malla fina reporta `Failed 3 mesh checks`:

| Check | Malla fina (cuarto) | Qué es |
|---|---|---|
| High aspect ratio (> 1000) | 3 030 celdas, máx 3026 | celdas del farfield (r > 2 m) justo detrás de la base: la primera celda axial de 0.5 mm sobre la base se extiende radialmente hasta el farfield, donde las celdas radiales miden 1.6 m. Flujo uniforme ahí; inofensivo. |
| Small determinant (< 0.001) | 529 k celdas | celdas estiradas (capa límite y farfield). Es un check de buena definición para solvers de movimiento de malla, no de `simpleFoam`. |
| Low interpolation weight (< 0.05) | 498 caras | idem, en los cambios de tamaño más fuertes |

Los checks que importan pasan: 100 % hexaedros, `Boundary openness OK`,
no-ortogonalidad máxima < 90° (86.7°, 7 563 caras > 70°), skewness 2.75 < 4.
`fvSchemes` ya lleva `limited corrected 0.33` y un corrector no-ortogonal por
esas 7 k caras.

## 3. Refinamiento local alrededor de las aletas (opcional)

Este paso no hace falta para la mayoría de las corridas. Antes de usarlo,
conviene saber **qué ya está refinado y qué no**:

| Dirección | Cómo se refina | ¿Es local? |
|---|---|---|
| cuerda (axial) | `FIN_H_X`, que parte `cyl` en `cyl` + `cylfin` | **sí**, sólo sobre la aleta |
| normal a la aleta (azimutal) | `AZ_FIN_H` y los bloques que tocan los planos | **sí**, sólo cerca de la aleta |
| envergadura (radial) | `FIN_H_R`, una zona radial | **no**: es un cilindro que va del inlet al outlet |

Las dos primeras ya son locales y no cuestan nada de más. La radial no puede
serlo dentro del generador: en una malla multibloque transfinita conforme, el
número de celdas de una dirección es el mismo a lo largo de toda la línea de
índice, así que una caja alrededor de la aleta necesitaría interfaces no
conformes. Por eso `FIN_H_R` sube el total un 36 % aunque sólo te interese
medio metro de dominio.

`./Allrefine` resuelve eso por afuera del generador, con `refineMesh` de
OpenFOAM: parte 2:1 las celdas de una región y nada más.

```bash
cd ~/runs/mi_corrida
./Allrefine            # el anillo de la punta (default)
./Allrefine tip        # lo mismo
./Allrefine fins       # la aleta entera, raíz a punta
```

o en un solo paso al crear la corrida:

```bash
./newCase.sh ~/runs/tip mesh/output/aconcagua_half_s3.msh --alpha 5 --refine tip
```

**Va entre `Allmesh` y `Allrun`, sobre una malla sin campos.** `refineMesh`
mapea los campos que encuentra, así que refinar después de correr te refina un
resultado a medias en lugar de la malla.

### Las dos regiones

Las dos son un cilindro alrededor del eje del cuerpo, y se diferencian sólo en
si llegan hasta la pared. La posición sale de la geometría de las aletas que
`build.py` escribe en `constant/meshInfo`, así que siguen a las aletas si
cambian. Los tres valores que la definen están arriba de
`system/topoSetDict` (`margin`, `wake`, `tipBand`).

| | `tip` | `fins` |
|---|---|---|
| región | anillo r de 0.186 a 0.260 m | el mismo cilindro sin el hueco |
| toca la capa límite | no | **sí** |

Medido sobre `coarse half`, 1 060 000 celdas:

| | `tip` | `fins` |
|---|---|---|
| celdas en el set | 24 784 (2.3 %) | 200 580 (18.9 %) |
| celdas resultantes | 1 233 488 (**+16 %**) | 2 464 060 (**+132 %**) |
| poliedros | 8 064 (0.65 %) | 7 392 (0.30 %) |
| no-ortogonalidad máx | 87.5°, sin cambio | 87.5°, sin cambio |
| skewness máx | 2.55, sin cambio | 2.55, sin cambio |
| y+ en la pared | sin cambio | **a la mitad bajo las aletas** |

Los patches se refinan solos: `fins` pasó de 10 984 a 13 192 caras con `tip`.

### Lo que cuesta

La malla **deja de ser 100 % hexaédrica**. La transición 2:1 no se puede
representar con hexaedros, así que aparecen poliedros. Son menos del 1 % y la
calidad medida no se mueve, pero es un cambio real respecto de lo que produce
el generador. Un detalle contraintuitivo: los poliedros viven en la
**superficie** del set, así que el anillo delgado genera casi tantos como una
región siete veces más grande.

`fins` llega hasta la pared del cuerpo y parte `y1` al medio ahí. El
`nutUSpaldingWallFunction` lo tolera, pero una descomposición de arrastre por
componente pasa a comparar patches a dos y+ distintos. Para eso, `tip`.

### Cuándo usarlo

Para arrastre a α = 0, probablemente no lo necesites: lo que manda es la
superficie de la aleta, y la cuerda y la normal ya están finas. Para el
vórtice de punta a ángulo de ataque, sí.

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

### El borde en escalera del patch `fins`

Mirando el patch en ParaView se ve un contorno escalonado en el borde de
ataque, el de fuga y la punta. Conviene saber qué es y qué no es, porque el
nombre "escalera" se asocia a un problema mucho peor del que hay acá.

**No es una escalera de pared.** En una malla cut-cell o en un `snappyHexMesh`
sin snap, lo que queda escalonado es la **superficie**, y entonces la
**normal** de la pared está mal: el esfuerzo de corte y la proyección de la
presión salen mal, y eso sí arruina la correlación. Acá la aleta es una placa
plana que vive en el plano, la deformación pone los nodos exactamente sobre
z = −t_half con residuo 0.0 nm, y toda cara del patch tiene la normal
correcta. Medido: el 28 % del área tiene la normal a más de 8° del plano, y
esa fracción **no cambia** al refinar (28.3 %, 27.6 %, 27.6 %), así que es el
bisel y el difuminado de la punta, geometría real, no un artefacto.

**Lo que sí está cuantizado es el contorno.** Una cara del plano de simetría
pasa a `fins` si **alguno** de sus cuatro nodos se movió, así que la aleta se
extiende hacia afuera hasta una celda en todo el perímetro. Ese es el escalón
que se ve.

Contra el área mojada exacta de la superficie, 376.0 cm² por semi-aleta:

| `FIN_H_X` / `FIN_H_R` | área de la malla | exceso | celdas |
|---|---|---|---|
| 5.0 / 12 mm (`coarse`) | 390.4 cm² | +3.8 % | 530 k |
| 2.5 / 6 mm | 383.7 cm² | +2.0 % | 993 k |
| 1.25 / 3 mm | 379.5 cm² | +0.9 % | 2.09 M |

Convergencia de primer orden limpia: al partir al medio los espaciamientos de
la aleta, el error se parte al medio.

**Qué le hace a los resultados.** El arrastre de fricción sobre las aletas es
proporcional al área mojada, así que se sobreestima un 3.8 % en `coarse`. Las
aletas son alrededor del 22 % del área mojada total, o sea menos de 1 % sobre
la fricción total y menos todavía sobre `Cd₀`: por debajo de la incertidumbre
del modelo de turbulencia. El arrastre de **presión** casi no se toca, porque
las caras del anillo quedan casi en el plano y su vector de área apunta a ±z,
perpendicular a la corriente. Lo que sí hereda el sesgo es `CN_α`, porque la
planform efectiva es un 3.8 % mayor, y con él un corrimiento chico del centro
de presión. Para correlacionar contra Barrowman eso es un sesgo conocido de
~4 % en `coarse` y ~1 % con la aleta al doble de resolución.

**Por qué la regla no puede ser "la cara es de la aleta si su centro lo es".**
Una cara con algún nodo sobre la aleta tiene ese nodo **duplicado** al coser
dos cuadrantes: son las dos caras de la aleta de espesor completo. Esa cara no
puede volverse interior. Etiquetarla `symm` pondría una condición de simetría
sobre lo que físicamente es la cuña que cierra la aleta, y dejaría un par de
caras de simetría casi coincidentes en medio del dominio. La regla del "algún
nodo" es lo que mantiene estancas las mallas de mitad y completa; el anillo de
una celda es su precio.

**Nota.** El preset `coarse` está hoy **mejor** que `fine` en la aleta
(+3.8 % contra +6.9 %), porque lleva `FIN_H_R = 12 mm`, que pone un límite de
zona justo en la punta y hace que la aleta cierre sobre una línea de nodos.
`fine` no lo trae. Si te importa ese sesgo en la malla de referencia,
`--preset fine --preset fintip`.

Convenciones del `forceCoeffs`:

- `dragDir` es la dirección de la corriente libre, `liftDir` es perpendicular
  en el plano de cabeceo x–y. `Cl` es la fuerza normal-al-viento en ese plano;
  la fuerza lateral (con `beta`) está en `forces1`.
- `Aref` es la **sección del cuerpo** (π D²/4 = 0.01791 m²) × la fracción de
  360° que hay en la malla (¼, ½, 1). Los coeficientes salen en su valor
  verdadero para cualquier sector; las **fuerzas** de `forces1` hay que
  multiplicarlas por 4, 2 o 1.
- `lRef` es el diámetro (0.151 m). Los momentos son respecto de la **punta**
  (`CofR (0 0 0)`), eje de cabeceo z. El centro de presión, medido desde la
  punta: `x_cp = -Cm · lRef / CN`, con `CN ≈ Cl` a ángulos chicos. Misma
  convención que Barrowman/OpenRocket.
- y+: el objetivo es que la mayoría de las paredes estén en 20 < y+ < 100.
  El `nutUSpaldingWallFunction` es continuo en y+ y tolera las colas
  (borde de ataque de las aletas, boattail). En las aletas el primer
  espaciamiento lo fija `AZ_FIN_H`, no `YPLUS_TARGET`.

## 7. Estudios típicos

**Cambiar de velocidad.** La malla se dimensiona para una velocidad: `U` en
`meshParams.py` fija `y1` para que el y+ dé en el objetivo, y correr a otra
velocidad deja el y+ fuera de rango sin avisar. Una malla por rango:

```bash
python build.py --preset medium --set U=272     # M 0.8
python build.py --preset medium --set U=612     # M 1.8
```

La tabla de `y1` contra velocidad está en
[PARAMETERS.md](PARAMETERS.md#flujo-la-malla-se-dimensiona-para-una-velocidad).
`Allrun` avisa si la `Uinf` del caso no se corresponde con la malla.

**Convergencia de malla.** `H_SCALE_WALL = False` mantiene `y1` (y por lo
tanto y+ y las wall functions) mientras todo lo demás se engrosa:

```bash
for s in 3 2 1.4 1; do python build.py --scale $s --set H_SCALE_WALL=False; done
```

**Barrido en ángulo de ataque.** Una malla, N corridas:

```bash
python build.py --preset medium --sector half
for a in 0 2 4 6 8; do
    ./newCase.sh ~/runs/sweep/a$a mesh/output/aconcagua_half_s2.msh --alpha $a --np 8
    (cd ~/runs/sweep/a$a && ./Allrun)
done
```

**Refinar las aletas.** Tres direcciones, tres parámetros (ver
[PARAMETERS.md](PARAMETERS.md#aletas)): `AZ_FIN_H` (normal a la aleta),
`FIN_H_X` (a lo largo de la cuerda), y `FIN_H_R` para la envergadura
(`--preset fintip`). La radial es global en x y cuesta +36 %; si sólo querés
la punta, `./Allrefine tip` hace lo mismo por +16 % (sección 3).

**Estela transitoria.** `--preset wake_unsteady --sector full`, y en el caso
cambiar a `pimpleFoam` (`ddtSchemes backward`, `deltaT` para Co ≈ 1 en la
estela cercana: con celdas de 6 mm y 100 m/s, ~5e-5 s). Un plano de simetría
suprime el desprendimiento asimétrico: usar la malla completa.

**Compresible.** A 100 m/s es M 0.29 e incompresible está bien. Para M > 0.3
la misma malla sirve con `rhoSimpleFoam`/`rhoPimpleFoam`: cambian `p` (Pa),
aparece `T` y `thermophysicalProperties`; la malla y `flowConditions` no.

## 8. Convenciones

- Los `.msh` no se versionan. Lo que se versiona es el `.params.py` (chico)
  junto a los resultados que produjo.
- Un directorio de corrida por condición, nombrado por lo que cambia
  (`a05_half_s2`, `ref_quarter_s1`). El `constant/meshInfo` de cada corrida
  dice de qué malla salió.
- `case/` no se corre. Si un cambio al caso vale para todas las corridas
  futuras, va a la plantilla y se commitea.
- Cambios de geometría van a `aconcaguaGeom.py` y se verifican con
  `python aconcaguaGeom.py` (ver [GEOMETRY.md](GEOMETRY.md)).

## Chuleta

```bash
# malla
python build.py --plan [--preset X] [--scale H] [--sector S]
python build.py --preset coarse
python build.py --preset medium --sector half
python build.py --preset fine --sector full --no-audit
python build.py --preset output/aconcagua_half_s2.params.py     # reconstruir

# caso
source /usr/lib/openfoam/openfoam2412/etc/bashrc
./newCase.sh ~/runs/NOMBRE mesh/output/MALLA.msh --alpha A --np N
./newCase.sh ~/runs/NOMBRE mesh/output/MALLA.msh --alpha A --refine tip
cd ~/runs/NOMBRE && ./Allrefine tip     # opcional, antes de Allrun
cd ~/runs/NOMBRE && ./Allrun
tail -n 3 postProcessing/forceCoeffs1/0/coefficient.dat
grep -A5 "yPlus yPlus write" log.simpleFoam | tail -5
```
