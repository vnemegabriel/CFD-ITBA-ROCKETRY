# Aconcagua — malla estructurada en gmsh: notas de diseño

*Por qué la malla está hecha así. Para usarla, [WORKFLOW.md](WORKFLOW.md); los
parámetros uno por uno, [PARAMETROS.md](PARAMETROS.md); medias y completas,
[SECTORES_Y_AOA.md](SECTORES_Y_AOA.md); la geometría, [GEOMETRÍA.md](GEOMETRÍA.md).*

Cuadrante 100 % hexaédrico, aletas resueltas, estela amplificada. Las medias
(2 cuadrantes) y completas (4) se arman rotando y cosiendo este cuadrante, con
la misma calidad por copia.

**Los números de referencia no están en este archivo.** Los genera
`python mesh/check.py --full` y viven en `mesh/check.baseline.json`, así que no
envejecen cuando algo cambia. Lo que sí está acá son las razones.

---

## 1. Zonas

Los dos ejes de la figura son lineales a trozos para que se vea cada bloque: no
está a escala.

![zonas meridionales](img/zones_meridional.png)

La línea cobre es el **shell 0**, la banda radial fina. Aguas abajo de la base
**se abre de 0.35 m a 0.70 m** para seguir a la estela que se ensancha, en vez
de que la estela se le salga. En ciruela, el núcleo butterfly y su anillo, que
existen sólo donde el eje es fluido: aguas arriba del casquete de nariz y aguas
abajo de la base.

![zonas en el plano transversal](img/zones_crossplane.png)

El cuarto de sección transversal lleva 12 bloques azimutales × 3 celdas, con los
anchos de bloque graduados hacia los dos planos de simetría porque ahí están las
aletas. El núcleo es una grilla de 18 × 18 sub-bloques: un núcleo cuadrado único
sólo puede presentarle dos aristas al anillo, así que `N_AZ_BLOCKS > 2` obliga a
un núcleo (N/2) × (N/2).

| Nivel | Zona | Hasta | Celdas radiales | Existe donde |
|---|---|---|---|---|
| L0 | núcleo | cuadrado de semiancho 0.45·Rᵢ | 18 × 18 | sólo con eje fluido |
| L1 | anillo | cuadrado → círculo interior Rᵢ(x) | 10 | sólo con eje fluido |
| L2 | zona de aleta | pared → `FIN_TIP_R` = 0.2355 m | 49 | siempre |
| L3 | zona 0 | → `ZONE_R[0]` = 0.35 m (0.70 m en el outlet) | 7 | siempre |
| L4 | zona 1 | → `ZONE_R[1]` = 2.00 m | 21 | siempre |
| L5 | zona 2 | → `ZONE_R[2]` = 11.775 m — **este es el farfield** | 14 | siempre |

A `H_SCALE 1`: y₁ = 303 µm, 31 celdas dentro de δ = 38 mm, 91 celdas radiales
en total, 144 azimutales alrededor.

| Segmento | x₀ | x₁ | Celdas axiales | Patch de pared |
|---|---|---|---|---|
| up (×6) | −17.730 | 0.026 | 80 | — |
| nose (×14) | 0.026 | 0.800 | 153 | cone |
| cyl | 0.800 | 2.029 | 102 | walls |
| cylfin | 2.029 | 2.530 | 63 | walls |
| finchord | 2.530001 | 2.830002 | 60 | walls |
| tail | 2.830002 | 2.955 | 25 | tail |
| wake1 | 2.955 | 3.255 | 110 | — |
| wake2 | 3.255 | 6.955 | 222 | — |
| wake3 | 6.955 | 41.370 | 91 | — |

`finchord` **es** la cuerda de la aleta: sus dos estaciones de borde están
dobladas sobre los bordes de ataque y de fuga, así que sus 60 celdas axiales
cubren 300 mm en la raíz y 150 mm en la punta — el índice axial *es* la fracción
de cuerda. `cylfin` es la aproximación, con la cara de aguas arriba plana y la
de aguas abajo sobre el borde de ataque.

---

## 2. Las aletas

![planform y sección de la aleta](img/zones_fin.png)

### Entran deformando, sin agregar bloques

La coordenada azimutal se deforma:

```
theta -> theta_f(x,r) + theta * (90 - 2 theta_f) / 90     theta_f = arcsin(t_half / r)
```

Un nodo que estaba en el plano de simetría aterriza en `z = -t_half`, que **es**
la superficie de la aleta: geometría exacta, topología intacta. Funciona sólo
porque la aleta real está biselada — el espesor va a cero de forma continua en
el borde de ataque y en el de fuga, así que la deformación se relaja sola ahí.
Los dos planos de simetría la reciben: el cuadrante contiene dos medias aletas,
una en cada plano. `theta_f` vale 2.28° en la raíz y 0.73° en la punta.

**Qué cara es aleta.** Una cara del plano de simetría con *algún* nodo desplazado
queda fuera del plano, así que es `fins`; sólo las que tienen los cuatro nodos
quietos siguen siendo `symm`. Esa regla no se puede hacer más lista: al ensamblar
sectores, una cara del plano se fusiona al interior únicamente si sus cuatro
nodos siguen en z = 0, así que "algún nodo desplazado → pared" está **forzado por
la topología**. Un test por centroide bajaría el sesgo de área de +10.3 % a
−0.2 % y rompería el cosido. La misma función decide las dos cosas, así que el
cuarto, la media y la completa no pueden discrepar (`sectorAssembly.py`).

Lo que hay que arreglar entonces no es la regla, es su entrada.

### Los bordes son frontera de bloque

Los bordes del planform tienen flecha: dx/dr = 1.5626 en el de ataque, 0.6249 en
el de fuga. Estaciones axiales que son **planos** los cruzan en diagonal, y la
clasificación cuantiza el contorno a la celda. El parche sale siendo el planform
**dilatado una celda**, nunca más chico:

| | área del parche vs planform | filas en envergadura | celdas de cuerda |
|---|---|---|---|
| H_SCALE 1 | +10.3 % | 50 | 81 |
| H_SCALE 3 | +18.9 % | 16 | 27 |
| H_SCALE 6 | +21.1 % | 8 | 14 |

Y como el espesor va a cero en el borde de ataque, ese anillo de caras de más
casi no tiene espesor: una pestaña serrada parada delante del borde real.

`FIN_EDGE_FIT` corta las estaciones axiales en el plano (x, r) — `x -> x + d(x, r)`,
aplicado al array de nodos terminado igual que la deformación de la aleta, sin
agregar un bloque ni snapear nada:

| Estación ancla | Corrimiento que lleva |
|---|---|
| `fin_x_start()` = 2.0292 | 0, fija |
| `FIN_LE_X_WALL` = 2.530001 | `x_LE(r) − x_LE(R_BODY)` → cae sobre el borde de ataque |
| `FIN_TE_X_WALL` = 2.830002 | `x_TE(r) − x_TE(R_BODY)` → cae sobre el de fuga |
| `X_BASE` = 2.9550 | `FIN_EDGE_TAIL_RELIEF ×` lo anterior |
| `X_BASE + X_WAKE_1` = 3.2550 | 0, fija |

Con eso un nodo **sobre** el borde de ataque tiene t = 0 exacto, nunca se
desplaza, y la regla `any` cae sola sobre el contorno verdadero: no queda nada
que cuantizar.

El borde de fuga no necesita estación propia porque ya la tiene: la raíz del TE
está en 2.830002 y la junta cilindro/boattail en 2.830000, a 2.5 µm. `cyl_end()`
mueve la junta esos 2.5 µm en vez de agregar una segunda estación. **No es
pedantería**: dejándola en 2.830000 cada nodo de esa estación conserva 0.7 µm de
semiespesor en lugar de cero, `any` se lleva la primera columna entera de caras
de `tail` al parche — 76 caras, +1.8 % de área — y el escalón se convierte en una
tira. La raíz del LE cae en 2.530001 = X_BASE − 0.425. Ninguno de los dos números
es casualidad, la aleta fue dibujada contra el cuerpo, pero `validate()` lo
afirma en vez de confiar, porque un cambio de `L_CYL` lo rompería en silencio.

**La pared no se mueve**: cada corrimiento se mide desde su propio valor en
r = R_BODY y `G.fin_edge_x()` congela por debajo, así que `d` se anula idéntica
sobre la superficie del cuerpo. `snap()` corre antes y su resultado en nanómetros
sobrevive.

**Lo que cuesta.** Una cara apoyada en el borde de ataque queda a 57.4° de la
dirección axial, porque eso *es* la flecha. No baja mientras las líneas radiales
sean círculos: el único camino es un collar envolviendo el borde del planform, o
sea otra topología. Medido a `H_SCALE 3`, cuadrante, `FIN_H_R` = 12 mm:

| | apagado | encendido |
|---|---|---|
| celdas | 371,920 | 371,920 |
| caras de aleta | 3,962 | 4,560 |
| área mojada por semialeta | +4.0 % | **−0.8 %** |
| caras fuera del planform real | ~1 celda alrededor | **0** |
| no-ortogonalidad máx | 74.92° | 74.92° |
| caras > 70° | 658 | 658 |
| caras > 40° | 3,999 | 65,947 |
| skewness máx | 4.578 | 4.578 |

El máximo no se mueve porque vive en el butterfly cap, no en la aleta. Lo que sí
crece es la banda entre 40 y 70: ese es el corte, y es el precio del borde.

El −0.8 % es la referencia, no la malla: `wetted_area_analytic()` integra desde
`FIN_ROOT_R` = 0.075, que está 0.5 mm adentro de la pared donde arranca la malla
(−0.40 %), y pone la superficie en y = r en vez del cilindro y = √(r² − t²) donde
la deformación realmente la envuelve (−0.24 %).

`Curve In Surface` de gmsh no es alternativa: las entidades embebidas sólo andan
con los algoritmos no estructurados Delaunay y HXT, así que una curva embebida y
un bloque transfinito se excluyen. El borde tiene que ser frontera de bloque.

### Lo que la flecha le cuesta al bloque de aproximación

Doblar la cara de aguas abajo sobre el borde de ataque **estira** ese bloque,
porque la de aguas arriba sigue siendo un plano y la cantidad de celdas es fija
a lo largo del índice radial:

```
estiramiento = 1 + (FIN_TIP_LE − FIN_ROOT_LE) / FIN_X_LEAD = 1 + 0.2508 / L
```

Con los 60 mm que `FIN_X_LEAD` tenía cuando el bloque era una banda común, contra
una flecha de 250 mm, eso da 5.1: 61 mm de largo en la raíz y 311 mm en la punta,
con 12 celdas uniformes que pasan de 5 mm a 26 mm. Los tamaños axiales en el
radio de la punta quedaban

```
5.0 mm (cyl)  →  25.9 mm (cylfin)  →  2.5 mm (finchord)
```

un escalón 5:1 para arriba y 10:1 para abajo, con la isla gruesa justo donde va
el choque del borde de ataque. **`H_SCALE` no lo arregla**: todos los tamaños
escalan juntos, así que la relación es invariante.

Tres cambios:

* `FIN_X_LEAD` es 500 mm, lo que deja el estiramiento en 1.50 y la celda en el LE
  de la punta en 7.7 mm.
* `cylfin` está **graduado** desde el tamaño del cilindro hasta `FIN_H_X`, y `cyl`
  ya no se fuerza a terminar en `FIN_H_X`. Uniforme sobre 500 mm serían 100
  celdas; graduado son 32, y caen contra el borde de ataque.
* `validate_params()` calcula `fin_lead_stretch()` y rechaza cualquier valor por
  encima de `FIN_LEAD_MAX_STRETCH` (2.0), nombrando el lead que necesitás.

Medido en `coarse`, cuadrante:

| | lead 60 mm | lead 500 mm graduado |
|---|---|---|
| celdas | 530,000 | **472,592** |
| celda axial en el LE de la punta | 25.9 mm | **7.7 mm** |
| no-ortogonalidad máx | 75.09° | 75.09° |
| caras > 70° | 1,165 | 1,165 |
| caras > 40° | 89,089 | 115,649 |
| media | 5.78° | 7.59° |
| skewness máx | 2.541 | 2.541 |

Menos celdas y 3.4× mejor espaciado en el LE, pagado en la banda de 40 a 70.

Efecto lateral: `cyl` ahora respeta el tamaño que le pide `SEGMENTS` (12 mm a
`H_SCALE` 1) en vez de graduarse en silencio hasta `FIN_H_X` a lo largo de todo
su largo. El cilindro medio quedó más grueso; ahí se fueron las 57 k celdas. Si
lo querés más fino, `SEGMENTS['cyl']['h_end']` existe para eso.

### La punta es un chaflán y no puede ser cuadrada

La aleta es una deformación azimutal. Una punta **cuadrada** necesita una cara a
r constante que cubra el espesor, con celdas del lado de afuera y ninguna del de
adentro en ese rango azimutal. Una estructura de bloques conforme no lo da:
adentro y afuera de `FIN_TIP_R` la parametrización azimutal difiere justo en
`theta_f`, así que la celda que cruza la punta queda cortada azimutalmente a lo
largo de su extensión radial, y ese corte **es** el chaflán. Haría falta un corte
topológico, como el borde de ataque necesitó ser frontera de bloque. No hay
perilla.

Lo que sí estaba mal y está arreglado: el chaflán corría **hacia afuera** de
`FIN_TIP_R`, así que la aleta llegaba a r = 0.2395 en vez de 0.2355 — 4 mm de
envergadura y 1.7 % de planform que el vehículo no tiene, justo donde nace el
vórtice de punta. Ahora cierra hacia adentro y `fin_zone_r()` pone la línea de
nodos **en** `FIN_TIP_R`, así que la aleta termina sobre un nodo. Sin esa línea
el contorno se vuelve a cuantizar: `coarse` con `FIN_H_R = None` da +4.1 % de
área mojada contra −0.9 % con ella.

El chaflán que queda es `max(FIN_TIP_SMEAR, celda radial en la punta)`, o sea que
un smear por debajo de `FIN_H_R` no compra nada. Medido en `coarse`:

| `FIN_H_R` | celda radial en la punta | chaflán | celdas/cuadrante | |
|---|---|---|---|---|
| `None` | 6.6 mm | 6.6 mm | 367,776 | sin nodo en la punta, contorno cuantizado |
| **0.012** | 12.1 mm | 12.1 mm | 459,392 | default, y `coarse` |
| 0.008 | 8.0 mm | 8.0 mm | 583,728 | |
| 0.006 | 6.0 mm | 6.0 mm | 688,432 | `fintip` |
| 0.004 | 4.0 mm | 4.0 mm | 884,752 | |
| 0.003 (smear 2 mm) | 3.0 mm | 3.0 mm | 1,048,352 | |
| 0.002 (smear 2 mm) | 2.0 mm | 2.0 mm | 1,329,744 | |

Ojo con lo que hace `FIN_H_R = 0.012` en `coarse`: la pila natural ahí ya daba
6.6 mm, así que la zona **engruesa** la punta a 12 mm y encima cuesta 92 k
celdas. Lo que compra es la línea de nodos, que vale más. Si te importa el
vórtice de punta, 0.006 te da las dos cosas.

### El refinamiento de volumen es otro trabajo

Nada de esto sirve para hacer las **celdas** alrededor de la aleta más chicas;
sólo las hace regulares. El refinamiento local 2:1 es `case-*/Allrefine`, que
corre `topoSet` + `refineMesh` sobre una región que `system/topoSetDict` lee de
`constant/meshInfo`. Los espaciados de cuerda (`FIN_H_X`) y normal a la aleta
(`AZ_FIN_H`) ya son locales; el radial es un cilindro que recorre todo el
dominio, y ése es el que `Allrefine` existe para arreglar.

**Coarsening y la aleta.** `H_SCALE_FIN = False` (presets `coarse` y `medium`)
exime a la aleta de `H_SCALE`: `FIN_H_X`, `AZ_FIN_H`, `FIN_H_R` y la cantidad de
celdas de los dos bloques azimutales que la tocan. Con todo escalado por 3 la
aleta de 12 mm quedaba de una celda de espesor con los bordes a dos celdas: eso
no es una aleta, es un bulto con su planform.

---

## 3. Parámetros — `meshParams.py`

Vos das **tamaños de celda**; las cantidades salen de `c = (L−h₁)/(L−h_N)` y
`n = 1 + ln(h_N/h₁)/ln c`. Las dos en forma cerrada, nada iterado.

`meshParams.py` tiene los defaults, y todo se puede pisar desde `build.py`
(`--preset`, `--scale`, `--sector`, `--set CLAVE=VALOR`) sin editar el archivo.
El conjunto de parámetros que se usó de verdad se escribe al lado de cada malla
como `<nombre>.params.py`, y ese archivo se puede volver a pasar con `--preset`.

### Grosería global — un solo número

```python
H_SCALE      = 1.0     # multiplica cada tamaño de celda, divide cada cantidad
H_SCALE_WALL = True    # False mantiene y1 (y por lo tanto y+) mientras el resto escala
```

Todo se vuelve a resolver desde los tamaños escalados, así que las
distribuciones siguen siendo correctas en vez de quedar diezmadas.
`build.py --scale <H>` lo pone desde la línea de comandos.

Para un estudio de convergencia de malla, `H_SCALE_WALL = False`: y+ se queda en
32 mientras el resto se engruesa.

### Zonas de refinamiento — el bloque que se edita

```python
ZONE_R = [0.35, 2.00, 11.775]    # donde TERMINA cada zona; LA ÚLTIMA ES EL FARFIELD
ZONE_H = [0.020, 0.200, 1.600]   # tamaño de celda en el borde EXTERIOR de cada zona
WAKE_ZONE_K   = 0                # qué zona se abre hacia la estela
ZONE0_R_WAKE  = 0.70             # radio exterior de esa zona en el outlet
WAKE_SPREAD_P = 0.5              # una estela turbulenta se abre como x^(1/2)
F_INLET    = 0.63                # tubo de blend en el inlet, como FRACCIÓN de la zona 0
F_WAKE_OUT = 0.50                # tubo de blend en el outlet, igual
```

El radio de farfield se declara una sola vez, como `ZONE_R[-1]`: los shells se
**derivan** de estas dos listas. Los tubos de blend son fracciones y no metros
para que no puedan crecer más que la zona que los contiene.

`validate_params()` corre en cada build y nombra lo que está mal:

| Regla | Error si se rompe |
|---|---|
| `ZONE_R` estrictamente creciente | `ZONE_R must increase: [...]` |
| `ZONE_R[0] > R_BODY` | la zona 0 estaría adentro del cuerpo |
| `ZONE_R[k] ≤ ZONE0_R_WAKE < ZONE_R[k+1]`, k = `WAKE_ZONE_K` | esa zona se cerraría, o se comería a la siguiente |
| `FIN_X_LEAD` deja el bloque de aletas sobre el cilindro | nombra la x resultante y el rango válido |
| `FIN_X_LEAD` contra la flecha | nombra el estiramiento, el tope y el lead que necesitás |
| la raíz del TE sobre la junta cilindro/boattail | `FIN_EDGE_FIT` no puede anclar el borde de fuga si no |
| `ZONE_H[k]` más chico que la zona k | `ZONE_H[k] is not smaller than zone k` |
| `0 < F_INLET, F_WAKE_OUT < 1` | fuera de rango |
| `CORE_FRAC < 1/√2` | el cuadrado del núcleo atravesaría su propio anillo |
| `N_AZ_BLOCKS` par y ≥ 2 | la grilla del núcleo queda indefinida |

### El resto

| Parámetro | Valor | Qué hace |
|---|---|---|
| `N_AZ_BLOCKS` | 12 | bloques azimutales por cuadrante, **tiene que ser par** |
| `N_AZ_CELLS` | 3 | celdas por bloque → 36 por cuadrante, 144 alrededor |
| `AZ_BLOCK_GROWTH` | 1.50 | relación de anchos, desde los planos de simetría hacia adentro |
| `AZ_FIN_H` | 5e-4 m | primera celda azimutal en r = R_BODY en el bloque que toca cada plano. `None` = uniforme |
| `AZ_RELAX` | `True` | relaja las dos de arriba hacia uniforme lejos de las aletas — ver § *La distribución azimutal se relaja lejos de las aletas* |
| `AZ_RELAX_X0` / `X1` / `R` | 0.80 / 2.029 / 0.35 m | dónde empieza y termina esa relajación, y hasta qué radio es completa |
| `N_RING` | 10 | celdas a lo ancho del anillo butterfly |
| `YPLUS_TARGET` | 32 | y⁺ en el **centro** de la primera celda; y₁ = 303 µm, 31 celdas dentro de δ |
| `SEGMENTS` | `h_start, h_end` | tamaño axial por segmento. `None` = continuar desde el vecino |
| `F_UP_INLET` | 0.50 | celda axial en el inlet, como fracción del primer sub-bloque (con tope `F_UP_MAX`) |
| `H_WAKE_BASE` | 5e-4 m | primera celda contra la base plana — es una pared |
| `CAP_R_FRAC` | 0.10 | radio del casquete butterfly / R_BODY. 7.55 mm en x = 26.3 mm |
| `CORE_FRAC` | 0.45 | semiancho del núcleo / Rᵢ. Tiene que quedar < 0.707 |
| `X_WAKE_1` / `X_WAKE_2` | 0.30 / 4.00 m | la estela cercana llega a 27 diámetros |
| `FIN_H_X` / `FIN_X_LEAD` | 0.005 / 0.500 m | celda axial sobre la cuerda, y la aproximación delante |
| `FIN_LEAD_MAX_STRETCH` | 2.0 | tope del estiramiento del bloque de aproximación |
| `FIN_H_R` | 0.012 m | celda radial en la punta; también fija el chaflán |
| `FIN_TIP_SMEAR` | 0.004 m | altura del chaflán de punta, hacia adentro de `FIN_TIP_R` |
| `FIN_SECTION` | `'wedge'` | `wedge` / `diamond` / `biconvex` / `naca` / `naca_te` |
| `UPSTREAM_L` / `DOWNSTREAM_L` | 6 / 13 | dominio, en largos de cuerpo |

Las constantes de geometría viven en `aconcaguaGeom.py`: cinco números definen
el cuerpo (`L_NOSE`, `L_CYL`, `L_TAIL`, `R_BODY`, `R_BASE`), todo lo demás se
deriva, y `validate()` afirma los invariantes.

---

## 4. Scripts

```
python mesh/build.py --plan                 # imprime el plan sin construir nada
python mesh/build.py --preset fine          # build -> snap -> aletas -> auditoría -> escritura
python mesh/build.py --preset smoke         # prueba de pipeline, ~25 s
python mesh/build.py --preset medium --sector full
python mesh/check.py [--full]               # la batería de regresión
cd <run>; ./Allmesh <archivo.msh>           # gmshToFoam -> meshInfo -> tipos -> checkMesh -> renumberMesh
```

| Archivo | Qué es |
|---|---|
| `aconcaguaGeom.py` | geometría de registro, reemplaza al STL. Se auto-verifica |
| `meshParams.py` | defaults, validación, el solver de cantidades. El único que normalmente editás |
| `presets/*.py` | conjuntos de parámetros con nombre, apilados sobre los defaults |
| `build.py` | la línea de comandos: presets, overrides, sector, salida, sidecars |
| `check.py` | la batería de regresión, contra `check.baseline.json` |
| `blockTools.py` | capa de bloques estructurados con memoización sobre gmsh |
| `buildHexBody.py` | la topología de bloques |
| `meshIO.py` | modelo gmsh → arrays numpy; escritor .msh v2.2; sidecars |
| `meshFinish.py` | el pipeline del cuadrante: construir, mallar, extraer, snap, ajustar bordes, deformar, clasificar |
| `finEdge.py` | el corte que pone los bordes de la aleta sobre fronteras de bloque |
| `finPatch.py` | deformación de aleta y clasificación de caras, sobre arrays |
| `sectorAssembly.py` | cuadrante → media / completa, rotando y cosiendo |
| `meshQuality.py` | las métricas de OpenFOAM, calculadas antes que OpenFOAM |
| `../case-*/fixPatchTypes.py` | **se corre justo después de gmshToFoam, no es opcional** |

Todo lo que pasa después de `gmsh.model.mesh.generate(3)` trabaja sobre arrays:
ni una llamada a la API de gmsh por nodo, una pasada vectorizada para el snap y
otra para la deformación, y un escritor que emite exactamente lo que gmshToFoam
consume. Sólo se escriben los nodos que referencia algún hexaedro, así que los
puntos de control de spline nunca llegan al archivo.

```
pip install gmsh numpy scipy
sudo apt install libglu1-mesa    # Linux/WSL: sin esto, import gmsh tira OSError: libGLU.so.1
```

El `.msh` del cuadrante `fine` pesa cientos de MB en ASCII: **generalo local, no
lo copies**, y usá `--out` para sacarlo de OneDrive. El sidecar `.params.py` lo
reconstruye.

`gmshToFoam` deja **todos** los patches como `type patch`. `symm` falla ruidoso;
`cone`, `walls`, `tail` y `fins` fallan **en silencio** — wall functions, yPlus y
forceCoeffs quedan mal sin decir nada. `fixPatchTypes.py` los pone desde la
entrada `patchTypes` de `meshInfo` y da error con cualquier nombre que no
reconozca.

Patches: `inlet` `outlet` `symm` `box` `cone` `walls` `tail` `fins`. La cantidad
de caras en patches tiene que dar exactamente igual a la de caras de borde —
`build.py` lo chequea y se niega a escribir si no — así que nada cae en
`defaultFaces`.

`Aref` referencia **la sección transversal del cuerpo**, no la semi-envergadura
de la aleta, por la fracción de 360° presente en la malla: 0.0044770 m² para el
cuadrante, `lRef` = 0.1510 m. Los dos, y la lista de patches de pared, los
escribe `build.py` en el sidecar `.meshInfo` que incluye `controlDict`:
`forceCoeffs` aborta con un nombre de patch que no encuentra, y una lista
mantenida a mano ya se desincronizó de la malla más de una vez.

---

## 5. Calidad

### Contra el snappyHexMesh que reemplazó

| | snappyHexMesh | esta malla (cuadrante) |
|---|---|---|
| hexaedros | 96.9 % | **100 %** |
| poliedros / prismas | 67,047 / 10,045 | **0 / 0** |
| celdas cóncavas | 6,026 | **0** |
| caras ilegales | 11 | **0** |
| celdas con volumen ≤ 0 | — | **0** |
| celdas dentro de δ | 12–20 | **31** |
| residuo en la pared | — | **0.00 nm** |
| no-ortogonalidad media | 6.81° | **5.35°** |
| máx aspect ratio | 24.7 | 3,026 (3,030 celdas en el farfield; 63 en el resto) |

### Auditoría vs checkMesh

`meshQuality.py` usa las definiciones de skewness de OpenFOAM (caras internas
normalizadas por la extensión de la cara en la dirección del sesgo, caras de
borde contra la proyección normal del owner). La no-ortogonalidad sí difiere,
porque acá los centros de celda son promedios de vértices y en checkMesh son
centroides pesados por volumen: ~76° acá contra ~87° en checkMesh sobre la misma
malla. **Los números de checkMesh son los que cuentan.**

### Lo que queda abierto

Con `-allGeometry`, checkMesh marca celdas de aspect ratio alto, todas en el
farfield (r > 2 m) justo detrás de la base, donde la primera celda axial de
0.5 mm contra la base se encuentra con celdas radiales de 1.6 m: un bloque
transfinito lleva una sola distribución axial a todos los radios. Es inofensivo
ahí, y sacarlo pediría una distribución axial que varíe con el radio, o sea otra
topología. `fvSchemes` lleva `limited corrected 0.33` y
`nNonOrthogonalCorrectors 1`.

**Un defecto diagnosticado y sin arreglar todavía.** Se ve en el render y está
medido; el trabajo no está hecho.

*El tamaño de celda se invierte al cruzar r = FIN_TIP_R sobre la nariz.*
Sólo el shell 0 se re-ajusta por estación (`buildHexBody._counts`); los shells
de afuera conservan su ratio global. El borde interior del shell 0 es la pared y
el exterior un cilindro fijo, así que hacia la nariz el span crece de 160 a
228 mm con las mismas 49 celdas y la misma primera celda de 303 µm:

| x [m] | r pared | span sh0 | ratio sh0 | última celda sh0 | 1ª celda sh1 | salto |
|---|---|---|---|---|---|---|
| 0.0263 (cap) | 0.0075 | 0.2279 | 1.0902 | 19.13 mm | 12.00 mm | **1.59** |
| 0.10 | 0.0203 | 0.2152 | 1.0885 | 17.78 mm | 12.00 mm | 1.48 |
| 0.40 | 0.0534 | 0.1821 | 1.0836 | 14.34 mm | 12.00 mm | 1.19 |
| 0.80 (cilindro) | 0.0755 | 0.1600 | 1.0798 | 12.11 mm | 12.00 mm | 1.01 |

No es que el ratio cambie mucho, es que **se invierte el signo**: las celdas
crecen, se achican 1.59× y vuelven a crecer, justo donde se para el choque de
proa a M ≥ 1.2.

Lo correcto es que el borde del shell 0 sea un **offset constante desde la
pared** (160 mm, que en el cilindro da exactamente 0.2355 y conserva la línea de
nodos de la punta de aleta). `zone_r_out(k, x)` ya es función de x, así que la
maquinaria está. La trampa: moviendo sólo ese borde, el span del shell 1 cambia
y el salto se muda a su borde exterior. Para que no cascadee tienen que seguir
la pared **todos** los bordes de zona, y ahí el farfield deja de ser un cilindro
— lo que toca el patch `box`, el blend del inlet y la apertura de la estela.

### La distribución azimutal se relaja lejos de las aletas

`AZ_BLOCK_GROWTH` y `AZ_FIN_H` agrupan celdas contra los dos planos de simetría
**porque ahí están las aletas**, y las aletas ocupan x ∈ [2.53, 2.93], el
13.5 % del cuerpo:

```
anchos de bloque [deg]: 2.17 3.25 4.87 7.31 10.96 16.44 16.44 10.96 7.31 4.87 3.25 2.17
```

Aplicada en toda estación, esa distribución daba en la nariz a x = 0.10 un arco
de celda de 0.134 mm contra 1.94 mm en el medio (max/min = 14.4) conviviendo con
una celda radial de 303 µm: relación de aspecto pésima donde a α = 0 el flujo es
axisimétrico y la resolución azimutal no aporta nada. No costaba celdas — son 36
por cuadrante igual — era cómo estaban repartidas.

Ahora `th` es función de la estación **y del radio**:

```
th_j(x, r) = th_fin_j + lam(x, r) (th_uni_j - th_fin_j)
lam(x, r)  = (1 - A(x)) min(1, AZ_RELAX_R / r)
```

`A(x)` es un smoothstep de 0 en `AZ_RELAX_X0` (fin de la nariz, x = 0.80) a 1 en
`AZ_RELAX_X1` (`fin_x_start()`, x = 2.029), así que **la banda de la aleta y todo
lo de aguas abajo conservan exactamente la malla que tenían**. Sobre la pared de
la nariz el arco de celda pasa de 0.134 / 1.94 mm a 0.885 mm parejo: max/min
14.4 → 1.00.

**El 1/r no es cosmético, es lo que hace que esto se pueda pagar.** Relajar sólo
en x — la lectura obvia — mueve todos los nodos que están en `th_j`, y en el
farfield ese movimiento es r·dth = 11.775 × 12.4° = **2.55 m de arco**, que hay
que soltar en los 1.23 m de cilindro disponibles: 64° de no-ortogonalidad de
pendiente media, 72° en el pico del smoothstep, en todo el campo exterior. No es
cuestión de ajustar nada; no hay largo de cilindro que alcance.

Con el 1/r el nodo se corre `AZ_RELAX_R·dth`, **el mismo arco a cualquier
radio**. Salen dos cosas de ahí. El corte axial deja de depender de r: 76 mm en
1.23 m, 5.3° en el pico. Y afuera de `AZ_RELAX_R` el patrón queda *trasladado*,
no abierto en abanico, que es la forma barata de moverlo: lo que se paga es la
inclinación de las líneas radiales, atan(`AZ_RELAX_R`·dth/r) — 12.2° apenas
afuera de `AZ_RELAX_R`, 4.3° en r = 1 m, 0.4° en el farfield. Adentro de
`AZ_RELAX_R` no hay inclinación ninguna, porque ahí la relajación es una
rotación pura.

Por eso `AZ_RELAX_R` es el borde exterior de la banda radial fina (`ZONE_R[0]`,
0.35 m) y no el radio del cuerpo: deja toda la capa límite y toda la envergadura
de la aleta del lado sin inclinar, y manda la banda de 12° al campo grueso.
Anclado en `R_BODY` la no-ortogonalidad media da 7.97°; en `ZONE_R[0]`, 7.65°,
contra 7.56° sin relajar (`coarse`).

Lo que se movió en el baseline `coarse` (la cuenta de celdas no se mueve: son
las mismas 36 por cuadrante, repartidas distinto):

| Métrica | Antes | Ahora | De dónde sale |
|---|---|---|---|
| arco de celda en la pared de la nariz, max/min | 43.34 | **3.00** | el arreglo. El 3.00 que queda es `H_SCALE_FIN = False`, que le deja 3 celdas a los bloques de la aleta y 1 al resto; a `H_SCALE 1` da 14.4 → **1.00** |
| no-ortogonalidad media | 7.56° | 7.65° | la inclinación radial afuera de 0.35 m |
| no-ortogonalidad máxima | 75.09° | 76.35° | casquete butterfly, ver abajo |
| skewness máxima | 2.541 | 2.176 | |
| volumen total | 6390.107 | 6390.663 | +87 ppm, ver abajo |

Los dos que subieron:

- **La máxima está en el casquete butterfly**, donde ya estaba. Ahí `lam = 1`,
  o sea que el anillo que pega el cuadrado del núcleo contra el círculo interior
  ve ahora arcos parejos en vez de agrupados, y ese mapeo cuadrado→círculo
  encaja 1.26° peor (2.2 % en `fine`, y ahí también 5.7 % más caras arriba de
  70°). Es el mismo bloque que ya daba el máximo de la malla; no aparece un
  lugar malo nuevo. Lo que lo arreglaría de verdad no es esto: es que la arista
  del núcleo se distribuya como `atan(v/a)` en vez de pareja, que es la
  distribución angular que el anillo le pide. Es otro trabajo — toca el
  butterfly, no la azimutal.
- **El volumen** sube 87 ppm porque el farfield es un **polígono inscripto** en
  el cilindro, y emparejar los arcos ahí (`lam ≈ 0.03` a r = 11.775) agranda un
  poco el polígono. Los nodos siguen sobre el cilindro con 0.032 mm de error.

`AZ_RELAX = False` reconstruye la malla anterior **exacta**, no parecida: con
`lam = 0` en todas partes `az_angles()` y `az_coefs()` devuelven lo de antes y
`c_axial()` vuelve a las rectas, así que las métricas de `coarse` dan los
dígitos del baseline viejo. Eso es la configuración `norelax` de `check.py`, y
está en el baseline: si algún día se mueve, lo que se movió no es la
relajación.

### Calidad azimutal

Dos cosas la producen, y las dos hacen falta.

**Las curvas radiales se canonicalizan de adentro hacia afuera antes de
crearlas.** Una curva radial es la arista 0 de un bloque anular y la arista 2 de
su vecino, así que los dos la piden en direcciones opuestas y la primera llamada
que llega al memo fija en qué sentido corre la progresión. Sin canonicalizar,
una columna azimutal recibe la pila agrupada contra la pared y todas las demás
reciben su recíproca: primera celda radial de 11.5 mm contra 30.5 mm en la misma
estación, con las celdas achicándose hacia afuera en vez de crecer.

**Spans angulares chicos.** Un cuadrilátero transfinito interpola el radio lineal
en índice, así que el error escala con la sagita arco-cuerda del bloque. El
sector más ancho es de 16.4° (sagita 0.0103) contra los 45° (0.0761) de antes.

### Restricciones del toolchain de las que depende el código

Cada una obliga algo en `blockTools.py` o `buildHexBody.py`. Cambiar ese código
sin respetarlas produce una malla que se construye y está mal.

| Restricción | Qué obliga |
|---|---|
| una entidad memoizada se crea una sola vez, en la dirección del primero que la pide | canonicalizar la dirección antes de crear cualquier curva |
| gmsh espacia los puntos transfinitos de un spline por **longitud de arco, no por parámetro** | los puntos de control no pueden dictar posiciones de nodo; muestrear denso y pasar la progresión real |
| un blend potencia `ξ^q` con 1 < q < 2 tiene **curvatura no acotada** en ξ = 0 | el blend de aguas arriba se ajusta por tangente: `Rᵢ² = r_cap² + 2 r_cap s_cap t + λt²` |
| un cuadrilátero transfinito se despega de un meridiano **curvo** | la ojiva se parte en 14 bloques con estaciones que equidistribuyen \|r″\|^½; el error cae como 1/n² |
| cada Point del modelo gmsh carga un nodo de malla | `meshIO.extract()` se queda sólo con los nodos que referencia algún hexaedro |
| un grupo físico de gmsh es **por superficie**, pero una cara de bloque lleva aleta y simetría a la vez | el parche de aleta se parte a nivel de elemento en el `.msh` |
| `c ** n` desborda cuando n·ln c pasa ~709 | `stack_sum` está protegido en logaritmos |
| `gmshToFoam` tipa **todos** los patches como `patch` | `fixPatchTypes.py` tiene que listarlos todos y saltear el header `FoamFile` |
| las entidades embebidas de gmsh sólo andan con Delaunay y HXT | un borde no se puede embeber en un bloque transfinito: tiene que ser frontera de bloque |
