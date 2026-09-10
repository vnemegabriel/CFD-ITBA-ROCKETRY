# Validación

Ninguna de las tres plantillas está validada. Este documento dice **contra qué
compararlas y con qué cuentas**, para que el resultado sea verificable y no una
opinión. Las fórmulas de referencia están escritas con sus constantes para que
puedas evaluarlas a mano.

Orden recomendado: primero el subsónico, que es donde hay referencias
analíticas confiables, y recién después los compresibles.

## 1. Chequeos que no cuestan nada y hay que hacer siempre

| Chequeo | Qué debe dar | Qué significa si falla |
|---|---|---|
| `Cl` a α = 0 en malla `half` o `full` | ≈ 0 | la costura entre cuadrantes o la deformación de la aleta no son simétricas |
| `CmRoll` a α = 0 | ≈ 0 | idem |
| y+ por patch | grueso de la distribución en 20 a 100 | ajustar `nLayers`, ver [TROUBLESHOOTING.md](TROUBLESHOOTING.md#y) |
| balance de `Cd` cuarto contra mitad | igual dentro de 1 % | `Aref` mal escalado, o el sector no es equivalente |
| linealidad de `CN` entre α = 2° y 6° | pendiente constante | no convergió, o y+ fuera de rango en las aletas |

## 2. Subsónico: contra Barrowman

La referencia son las fórmulas de Barrowman, que son las que usan OpenRocket y
el resto de la cohetería amateur. Referidas a la **sección del cuerpo**
`A = πD²/4` y por radián.

- Barrowman (1967), *The Practical Calculation of the Aerodynamic
  Characteristics of Slender Finned Vehicles*, tesis de maestría, Catholic
  University of America. Es el documento original.
- Barrowman y Barrowman (1966), *A Method for Calculating the Static Margin of
  a Slender Missile*.
- Para la parte de aletas y su interferencia, la formulación equivalente está
  en Fleeman, *Tactical Missile Design*, 2ª ed., cap. 2.

### 2.1 Fuerza normal y centro de presión

Con `D = 0.151 m`, `L = 2.955 m`:

| Componente | `CN_α` [1/rad] | `x_cp` desde la punta [m] |
|---|---|---|
| ogiva | 2.0 (exacto para cualquier nariz de revolución) | 0.466 · L_nariz = 0.37 |
| cilindro | ≈ 0 en teoría esbelta | — |
| boattail | `2[(d_base/d)² − 1]` = −0.94 | 2.89 |
| aletas (4, con interferencia) | ver abajo | 2.70 |

El término de aletas de Barrowman, para `N = 4`:

```
CN_fin = (4 N (s/d)²) / (1 + sqrt(1 + (2 l / (c_r + c_t))²))
K_interf = 1 + r / (s + r)
```

con `s` = semi-envergadura expuesta = 0.2355 − 0.0755 = 0.160 m,
`r` = radio del cuerpo = 0.0755 m, `d` = 0.151 m,
`c_r` = 0.301 m, `c_t` = 0.150 m, `l` = longitud del borde de ataque.

**Cómo compararlo con la CFD.** Corré un barrido en α de 0 a 6° en malla
`half`, ajustá la pendiente de `Cl` contra α en radianes, y comparala con la
suma. Para el desglose por componente, agregá un `forceCoeffs` por patch
(`cone`, `walls`, `tail`, `fins`) y compará cada término por separado: es
mucho más informativo que un solo número.

Ojo con el sesgo conocido: el contorno del patch `fins` está cuantizado y su
área queda 3.8 % alta en `coarse`. Eso se traslada casi directo a `CN_fin`.
Ver [WORKFLOW.md](WORKFLOW.md#el-borde-en-escalera-del-patch-fins).

### 2.2 Arrastre a α = 0

Descomposición, que es lo que conviene comparar en vez del total:

| Término | Referencia |
|---|---|
| fricción | `Cf · S_mojada / Aref` con `Cf = 0.0576 Re_L^-0.2` (placa plana turbulenta, Schlichting *Boundary-Layer Theory* 7ª ed. cap. 21). A `Re_L = 2×10⁷` da `Cf ≈ 0.00203` |
| factor de forma del cuerpo | `1 + 1/(2·fineness)`, con fineness = L/D = 19.6 |
| factor de forma de la aleta | `1 + 2t/c` |
| arrastre de base | correlación de Hoerner, *Fluid-Dynamic Drag* (1965), cap. 3, con `(d_base/d)² = 0.53` |

La parte viscosa y la de presión salen separadas en `postProcessing/forces1`,
columna por columna, así que cada término se compara con su referencia.

## 3. Transónico y supersónico

Acá no hay fórmula analítica que valga como referencia. Tres caminos, de menos
a más esfuerzo:

**Contra el propio tutorial.** Antes de creerle al caso del cohete, corré el
tutorial del que sale la plantilla y reproducí su resultado. Es la forma más
barata de separar un problema de configuración de uno de malla:

```
$FOAM_TUTORIALS/compressible/rhoSimpleFoam/aerofoilNACA0012
$FOAM_TUTORIALS/compressible/rhoCentralFoam/biconic25-55Run35
```

El biconico trae datos experimentales, así que da una validación de verdad del
solver y del esquema de flujo.

**Contra soluciones analíticas de onda.** Para verificar que el esquema captura
bien un choque, la relación de Rankine-Hugoniot en una cuña es exacta:

```
tan(theta) = 2 cot(beta) (M² sin²(beta) − 1) / (M² (gamma + cos 2 beta) + 2)
```

Anderson, *Modern Compressible Flow*, 3ª ed., cap. 4. El tutorial
`rhoCentralFoam/wedge15Ma5` es justamente eso.

**Contra la curva `Cd(M)` de OpenRocket.** Es la comparación que más te importa
para la trayectoria, pero es la más débil como validación: OpenRocket usa
correlaciones semi-empíricas, no verdad de campo. Sirve para detectar un error
de un factor 2, no para justificar un 5 %.

### El pico transónico

Lo que hay que reproducir cualitativamente es la subida de `Cd` entre M 0.8 y
M 1.2, que para un cuerpo esbelto con aletas es del orden de un factor 2 a 3
respecto del valor subsónico. Si tu curva no lo tiene, o el `transonic yes` no
está haciendo efecto, o la malla es demasiado gruesa donde se forma el choque.

## 4. Convergencia de malla

Independiente del régimen, y necesaria para poder citar cualquier número:

Cambiá `maxCellSize` en `mesh/system/meshDict` por factores de 2 y dejá los
niveles quietos: eso escala todo el campo sin tocar la resolución relativa
entre zonas. Una malla por valor, el mismo punto de Mach en todas.

Ojo con y+: el espesor de la primera capa sigue a la celda de superficie, así
que se mueve con `maxCellSize`. Compensá con `nLayers` para que y+ quede en la
misma banda en las cuatro mallas, o el estudio mezcla dos efectos. Graficá `Cd`
contra `N^(-2/3)` y extrapolá; la pendiente te dice el orden observado y la
ordenada al origen el valor extrapolado. El método está en Roache (1994),
"Perspective: A Method for Uniform Reporting of Grid Refinement Studies",
*J. Fluids Eng.* **116**, 405-413.

## 5. Registro

Para cada corrida que quieras citar, guardá junto a los resultados:

- el `output/<nombre>.params.py` de la malla, que la reconstruye exacta
- el `constant/meshInfo` de la corrida
- el `system/flowConditions`
- el `log.checkMesh` y el y+ por patch del log del solver

Con eso cualquiera reproduce el número. Sin eso, no.
