# Parámetros

Referencia de todo lo que se puede cambiar sin tocar código: los parámetros
de la malla (`mesh/meshParams.py`, presets, `--set`), las condiciones del caso
(`case/system/flowConditions`) y lo que la malla le cuenta al caso
(`constant/meshInfo`).

## Malla — `meshParams.py`

Se dan **tamaños de celda**; los conteos se derivan del tramo que tienen que
llenar (`c = (L−h₁)/(L−h_N)`, `n = 1 + ln(h_N/h₁)/ln c`, en forma cerrada).
Nada se itera ni se "relaja" por calidad: "0.3 mm en la pared creciendo a
20 mm en r = 0.35 m" es una afirmación sobre la malla, no un deseo.

Toda variable de la tabla se puede fijar desde un preset o con
`--set NOMBRE=valor`. Un nombre que no está en `OVERRIDABLE` es un error.

### Grosor global

| Parámetro | Default | Qué hace |
|---|---|---|
| `H_SCALE` | 1.0 | multiplica todo tamaño de celda y divide todo conteo. Celdas ≈ 1/H³. `--scale` |
| `H_SCALE_WALL` | `True` | `False` mantiene fija **la primera celda de cada pared** mientras el resto escala: `y1` (y por lo tanto y+) en las paredes laterales, la celda axial en la punta del cono (`up.h_end`/`nose.h_start`) y la celda axial sobre la base (`H_WAKE_BASE`). Lo correcto para un estudio de convergencia con wall functions. La de la punta importa también por calidad: el casquete es casi vertical en el ápice y las celdas del núcleo miden ~0.1 mm, así que una celda axial de 4.5 mm ahí (`H_SCALE 3`) lleva el skewness de las caras de pared a 4.6 |
| `H_SCALE_FIN` | `True` | `False` mantiene la **aleta** a resolución completa mientras el resto escala: `FIN_H_X` (cuerda), `AZ_FIN_H` y las celdas de los dos bloques azimutales que tocan la aleta (normal), `FIN_H_R` (envergadura). Con todo escalado ×3 la aleta de 12 mm queda de una celda de espesor y sus bordes caen dos celdas aparte: no es una aleta, es un bulto. `coarse` y `medium` lo ponen en `False` |

### Flujo: la malla se dimensiona para una velocidad

| Parámetro | Default | Qué hace |
|---|---|---|
| `U`, `NU`, `RHO` | 100, 1.5e-5, 1.225 | fijan Re_L y con él `y1`. **No** son las condiciones del caso: esas van en `flowConditions` |
| `YPLUS_TARGET` | 32 | y+ en el **centro** de la primera celda → `y1 = 303 µm`; 31 celdas dentro de δ en la malla fina |

**La malla es específica de una velocidad, y esto es fácil de pasar por alto.**
`U` no entra en la simulación: solo sirve para resolver `y1` de manera que el
y+ dé en el objetivo. Si después corrés el caso a otra velocidad, la primera
celda sigue siendo la misma y el y+ se va de rango **sin ningún error**. El
sesgo es del tamaño del cambio de velocidad:

| `U` [m/s] | Mach | `y1` para y+ = 32 | y+ que da la malla de 100 m/s |
|---|---|---|---|
| 50 | 0.15 | 566 µm | 17 |
| **100** | **0.29** | **303 µm** | **32** |
| 170 | 0.50 | 188 µm | 52 |
| 272 | 0.80 | 123 µm | 79 |
| 340 | 1.00 | 101 µm | 96 |
| 612 | 1.80 | 59 µm | 163 |

La regla es: **una malla por rango de velocidad.** Se construye pasando la
velocidad de la corrida, y `y1` se re-resuelve solo:

```bash
python build.py --preset medium --set U=272     # para correr cerca de M 0.8
```

`Allrun` compara la `Uinf` de `flowConditions` con la `Uref` que la malla
dejó grabada en `constant/meshInfo`, estima el y+ resultante y avisa si se
fue de rango. Es un aviso, no un error: correr a otra velocidad es válido
mientras sepas qué y+ te queda.

### Zonas radiales (lo que más se edita)

Zonas concéntricas desde la pared hacia afuera. `ZONE_R[k]` es donde **termina**
la zona k, así que la lista crece y **la última entrada es el radio del
farfield**. `ZONE_H[k]` es el tamaño de celda en el borde **exterior** de la
zona k; el borde interior de la zona 0 es `y1`.

| Parámetro | Default | Qué hace |
|---|---|---|
| `ZONE_R` | `[0.35, 2.00, 11.775]` m | fin de cada zona; la última es el farfield |
| `ZONE_H` | `[0.020, 0.200, 1.600]` m | tamaño de celda al final de cada zona |
| `WAKE_ZONE_K` | 0 | qué zona se abre en la estela. Las zonas interiores a ella (la de `FIN_H_R`, por ejemplo) se abren en proporción, así ninguna queda aplastada contra el tubo de salida |
| `ZONE0_R_WAKE` | 0.70 m | radio exterior de esa zona en el outlet |
| `WAKE_SPREAD_P` | 0.5 | la zona se abre como x^p (una estela turbulenta ensancha como x^½) |
| `F_INLET` | 0.63 | radio del tubo interior en el inlet, como **fracción** de la zona más interior (la de `FIN_H_R` si existe, si no `ZONE_R[0]`) |
| `F_WAKE_OUT` | 0.50 | idem en el outlet, fracción del radio de la zona más interior ahí |

Agregar una zona es agregar una entrada a las dos listas. Ejemplo, una zona
de 6 mm que termine en r = 0.28 m:

```python
ZONE_R = [0.28, 0.35, 2.00, 11.775]
ZONE_H = [0.006, 0.020, 0.200, 1.600]
WAKE_ZONE_K = 1        # la zona que se abre en la estela ahora es la 1
```

Para la punta de las aletas hay un atajo: `FIN_H_R` (tabla de aletas) inserta
la zona en el radio exacto de la punta sin tocar las listas ni `WAKE_ZONE_K`.
Es lo que hace `presets/fintip.py`.

### Axial

| Parámetro | Default | Qué hace |
|---|---|---|
| `SEGMENTS` | ver archivo | `(nombre, dict(h_start, h_end))` por tramo: `up`, `nose`, `cyl`, `tail`, `wake1`, `wake2`, `wake3`. `None` = continuar del vecino, así no puede haber un salto en una interfaz de bloque |
| `F_UP_INLET`, `F_UP_MAX` | 0.50, 0.60 | celda axial en el inlet como fracción del primer sub-bloque aguas arriba (1.48 m en la fina) |
| `H_WAKE_BASE` | 5e-4 m | primera celda axial sobre la base plana (es pared) |
| `X_WAKE_1`, `X_WAKE_2` | 0.30, 4.00 m | fin de `wake1` y `wake2` detrás de la base (27 diámetros de estela cercana) |
| `UPSTREAM_L`, `DOWNSTREAM_L` | 6, 13 | dominio en longitudes de cuerpo: inlet en −17.7 m, outlet en +41.4 m |
| `N_NOSE_BLOCKS`, `N_UP_BLOCKS` | 14, 6 | sub-bloques del cono y del tubo aguas arriba, para mantener el error de interpolación transfinita bajo 2 % de `y1` |

### Azimutal

| Parámetro | Default | Qué hace |
|---|---|---|
| `N_AZ_BLOCKS` | 12 | bloques azimutales por cuadrante, **par** (el núcleo es una grilla (N/2)²) |
| `N_AZ_CELLS` | 3 | celdas por bloque → 36 por cuadrante, 144 en la vuelta |
| `AZ_BLOCK_GROWTH` | 1.50 | relación de anchos entre bloques sucesivos, desde los planos de simetría hacia 45° |
| `AZ_FIN_H` | 5e-4 m | primera celda azimutal en r = R_BODY en el bloque que toca cada plano: **es el espaciamiento normal a la aleta** (y+ ≈ 50 en la aleta). `None` = uniforme |
| `N_RING` | 10 | celdas a través del anillo del butterfly |
| `CAP_R_FRAC` | 0.10 | radio del casquete butterfly sobre el cono / R_BODY |
| `CORE_FRAC` | 0.45 | semiancho del núcleo cuadrado / radio interior local; < 0.707 |

### Aletas

| Parámetro | Default | Qué hace |
|---|---|---|
| `FINS_ON` | `True` | `--fins` / `--no-fins` |
| `FIN_SECTION` | `'wedge'` | `wedge` (como está dibujada: biseles absolutos), `diamond`, `biconvex`, `naca`, `naca_te` |
| `FIN_TIP_SMEAR` | 0.004 m | banda radial en la que el espesor cierra a cero en la punta |
| `FIN_H_X` | 0.005 m | celda axial sobre la cuerda; parte `cyl` en `cyl` + `cylfin` y refina `tail`. `None` = sin refinar |
| `FIN_X_LEAD` | 0.060 m | cilindro incluido delante del borde de ataque de la raíz en el bloque refinado |
| `FIN_H_R` | `None` | celda **radial en la punta** de la aleta. `None` = sin zona propia: la punta recibe lo que da la pila de zonas (12 mm en `fine`, 40 mm a `H_SCALE 3`). Un valor inserta un límite de zona en `FIN_TIP_R + FIN_TIP_SMEAR` (la aleta cierra sobre una línea de nodos) y la pila radial crece de `y1` a `FIN_H_R` sobre la envergadura. Se inserta en su radio dentro de `ZONE_R`; `WAKE_ZONE_K` sigue indexando **tu** lista. `coarse`/`medium` usan 0.012, `fintip` 0.006 |

La planform, el espesor y los biseles están en `aconcaguaGeom.py`
([GEOMETRY.md](GEOMETRY.md)): son geometría, no malla.

### Sector

| Parámetro | Default | Qué hace |
|---|---|---|
| `SECTOR` | `'quarter'` | `quarter`, `half`, `full`. `--sector`. Ver [SECTORS_AND_AOA.md](SECTORS_AND_AOA.md) |

### Validación

`validate_params()` corre antes de cualquier build y nombra lo que está mal:

| Regla | Mensaje |
|---|---|
| `ZONE_R` estrictamente creciente | `ZONE_R must increase` |
| `ZONE_R[0] > R_BODY` | `ZONE_R[0] = ... is inside the body` |
| `ZONE_R[k] ≤ ZONE0_R_WAKE < ZONE_R[k+1]` con k = `WAKE_ZONE_K` | la zona se cerraría, o se tragaría a la siguiente en el outlet |
| `ZONE_H[k]` menor que el ancho de la zona k | `ZONE_H[k] ... is not smaller than zone k` |
| `FIN_X_LEAD` deja el inicio del bloque de aletas sobre el cilindro | nombra la x resultante y el rango válido |
| `0 < F_INLET, F_WAKE_OUT < 1`; `H_SCALE > 0`; `CORE_FRAC < 1/√2`; `N_AZ_BLOCKS` par | fuera de rango |
| la zona de `FIN_H_R` no cae a menos de dos celdas de una entrada de `ZONE_R` | `FIN_H_R puts a zone boundary at r = ... within two cells of ZONE_R entry ...` |
| `SECTOR`, `FIN_SECTION` conocidos | |

### Derivados (no se editan)

`r_far()`, `shell_spec()`, `derived()` (Re, u_τ, y1, δ, celdas por shell),
`axial_plan()`, `az_angles()`, `predicted_cells()`. `python meshParams.py`
o `build.py --plan` los imprimen.

## Caso — `system/flowConditions`

Hay una versión por plantilla. La subsónica se parametriza por velocidad; las
dos compresibles por **número de Mach** más el estado ambiente, y derivan la
velocidad. Modelos y referencias en [SOLVERS.md](SOLVERS.md).

### Compresibles (`case-transonic/`, `case-supersonic/`)

| Entrada | Default | Qué es |
|---|---|---|
| `Minf` | 0.80 / 1.80 | número de Mach de la corriente libre |
| `pInf` | 101325 | Pa, presión estática ambiente (ISA nivel del mar) |
| `Tinf` | 288.15 | K, temperatura estática ambiente |
| `alpha`, `beta` | 0 | deg, igual que en la subsónica |
| `Ti`, `nuRatio` | 0.005, 5 | igual que en la subsónica |

`Uinf`, `rhoInf`, `muInf` y `nuInf` salen derivadas en `flowDerived`; las
fórmulas y su verificación numérica están en
[SOLVERS.md](SOLVERS.md#5-de-qué-se-derivan-las-condiciones).

### Subsónica (`case-subsonic/`)

| Entrada | Default | Qué es |
|---|---|---|
| `Uinf` | 100 | m/s |
| `alpha` | 0 | deg, ángulo de ataque; U rota de +x hacia +y. Requiere `half` o `full` |
| `beta` | 0 | deg, deslizamiento; hacia +z. Requiere `full` |
| `nu` | 1.5e-5 | m²/s |
| `rhoInf` | 1.225 | kg/m³, sólo para escalar coeficientes y fuerzas |
| `Ti` | 0.005 | intensidad de turbulencia en la corriente libre → `k = 1.5 (U·Ti)²` |
| `nuRatio` | 5 | `nut/nu` en la corriente libre → `omega = k / nut` |

`system/flowDerived` deriva `Ux Uy Uz`, `dx dy dz` (drag), `lx ly lz`
(lift), `kInlet`, `omegaInlet`, `nutInlet`. Los campos de `0.orig/` sólo
referencian esos nombres.

## Refinado local — `system/topoSetDict`

Sólo lo usa `./Allrefine`, el paso opcional que refina 2:1 alrededor de las
aletas. Los tres valores de arriba del archivo definen la región; el resto se
deriva de la geometría de las aletas que viaja en `constant/meshInfo`, así que
la región sigue a las aletas si cambian.

| Entrada | Default | Qué es |
|---|---|---|
| `margin` | 0.020 m | holgura alrededor de la planform, en x y en r |
| `wake` | 0.050 m | cuánto se extiende la región detrás del borde de fuga, para seguir el vórtice de punta |
| `tipBand` | 0.050 m | espesor del anillo hacia adentro desde la punta. Sólo lo usa `finTip` |

Los dos sets que salen de ahí:

| Set | Región | Toca la pared |
|---|---|---|
| `finTip` | anillo entre `finTipR − tipBand` y `finTipR + finTipSmear + margin` | no, y+ intacto |
| `finBox` | el mismo cilindro sin el hueco: aleta entera, raíz a punta | sí, parte `y1` al medio |

`system/refineMeshDict` casi no se toca. `set` lo sobreescribe `Allrefine`
según su argumento, y `directions` conviene dejarlo en las tres: el sistema de
coordenadas de `refineMesh` es **cartesiano** y la malla es cilíndrica, así que
refinar sólo dos direcciones refina un subconjunto que rota con el azimut.

## Lo que la malla le dice al caso — `constant/meshInfo`

Escrito por `build.py` junto al `.msh`; `Allmesh` lo copia al caso.

| Entrada | Ejemplo | Quién lo usa |
|---|---|---|
| `sector`, `symmetryFraction` | `half`, `0.5` | `Allrun` (chequeo de ángulos), vos |
| `Aref`, `lRef` | `0.00895 m²`, `0.151 m` | `forceCoeffs` |
| `wallPatches` | `(cone walls tail fins)` | `forceCoeffs`, `forces`, `wallShearStress` |
| `patchTypes` | `{symm symmetry; cone wall; ...}` | `fixPatchTypes.py` |
| `finTipR`, `finTipSmear`, `finX0`, `finX1`, `finRootR` | `0.2355`, `0.004`, `2.52922`, `2.93` | `topoSetDict` (`./Allrefine`) |
| `nCells`, `hScale`, `y1`, `yPlusTarget`, `fins` | | registro. `Allrefine` actualiza `nCells` |
