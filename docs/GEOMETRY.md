# Geometría

`mesh/aconcaguaGeom.py` es la **geometría de registro** del Aconcagua. No hay
STL ni STEP en el pipeline: el cuerpo es una fórmula y la malla se apoya en
ella directamente. Un STL es una muestra de 2 096 facetas de estas mismas
fórmulas (las constantes se midieron sobre él y se verificaron contra la
forma analítica con 5 cifras).

Por qué importa para mallar: una superficie triangulada no tiene plano
tangente ni curvatura, y todo mallador que la consume tiene que adivinarlos
(`surfaceFeatureExtract`, `resolveFeatureAngle`, tolerancias de snap...). Con
`r(x)` en forma cerrada se deriva en vez de adivinar: normales exactas, y un
nodo se **pone** sobre la superficie en vez de iterarse hacia ella.

## Los cinco números

```python
L_NOSE   = 0.80000   # m  longitud del cono (Von Kármán / LV-Haack, C = 0)
L_CYL    = 2.03000   # m  cilindro
L_TAIL   = 0.12500   # m  boattail (tronco de cono recto)
R_BODY   = 0.07550   # m  radio del cilindro  (D = 0.151)
R_BASE   = 0.05500   # m  radio de la base
```

Todo lo demás se deriva: `X_BODY_1 = L_NOSE`, `X_BODY_2 = X_BODY_1 + L_CYL`,
`X_BASE = L_TOTAL = 2.955`, `D_BODY`, `BOATTAIL_HALF_ANGLE`. Derivar las
uniones en vez de declararlas es lo que impide que un cambio en una longitud
abra un escalón en `r_body` sin que nadie lo note.

Perfil: `r_body(x)` (ogiva | cilindro | cono), `drdx_body(x)`,
`wall_normal(x)`, `x_nose_of_r(r)` (inversa de la ogiva, por bisección: la
pendiente en la punta es **infinita**, `r ~ x^¾`, y eso es lo que dicta la
topología de casquete butterfly sobre el cono).

## Aletas

Cuatro deltas recortadas planas a 0°, 90°, 180°, 270°, ancladas a la **base**
para que viajen con la cola si cambia el largo del cuerpo:

```python
N_FINS       = 4          # fijado por la topología de la malla, no se cambia
FIN_T        = 0.012      # m  espesor total
FIN_ROOT_R   = 0.07500    # m  radio de la cuerda raíz (enterrada 0.5 mm en el cuerpo)
FIN_TIP_R    = 0.23550    # m  semi-envergadura
FIN_ROOT_LE  = X_BASE - 0.42578
FIN_ROOT_TE  = X_BASE - 0.12531    # = X_BODY_2: la raíz termina donde empieza el boattail
FIN_TIP_LE   = X_BASE - 0.17500
FIN_TIP_TE   = X_BASE - 0.02500    # la punta termina sobre el boattail
FIN_LE_BEVEL = 0.03320    # m  bisel del borde de ataque (ABSOLUTO, no % de cuerda)
FIN_TE_BEVEL = 0.02116    # m  bisel del borde de fuga
```

`fin_planform(σ)` da borde de ataque y cuerda a la fracción de envergadura σ;
`fin_halfthickness(ξ, cuerda, section)` el semiespesor a la fracción de cuerda
ξ para las secciones `wedge` (como está dibujada), `diamond`, `biconvex`,
`naca`, `naca_te`; `fin_half_thickness(x, r)` es el campo que lee la
deformación de la malla. Los biseles son absolutos, así que la sección cambia
de la raíz (11 % de la cuerda) a la punta (22 %): es lo que hace una aleta
mecanizada.

La malla introduce las aletas **deformando la coordenada azimutal**, sin
agregar bloques, y eso funciona sólo porque el espesor va a cero de forma
continua en los bordes (bisel). Una sección con borde de ataque redondo
(`naca`) sigue siendo analítica pero la deformación en el borde ya no es
suave; para resolverla de verdad haría falta una topología C alrededor del
borde. Ver [MESH_DESIGN.md](MESH_DESIGN.md#how-the-fins-went-in).

## Cambiar la geometría

1. Editar las constantes en `aconcaguaGeom.py`. Nada más: las uniones, las
   estaciones de las aletas y los valores de referencia siguen solos.
2. `python aconcaguaGeom.py` corre `validate()`: cierre del cono sobre el
   cilindro, monotonía de la ogiva (el mapa del casquete debe ser biyectivo),
   raíz de la aleta dentro del cuerpo, biseles que no se solapan en la punta,
   `N_FINS == 4`. Si hay un STL, `python aconcaguaGeom.py ruta.stl` compara
   `r_body` contra sus vértices.
3. `python build.py --plan`: `validate_params()` de la malla vuelve a chequear
   los acoples con la geometría (`ZONE_R[0] > R_BODY`, `FIN_X_LEAD` sobre el
   cilindro, etc.) y nombra lo que hay que ajustar.

Qué parámetros de malla suelen moverse con la geometría:

| Cambio | Revisar |
|---|---|
| `R_BODY` mayor | `ZONE_R[0]` (tiene que quedar fuera del cuerpo), `ZONE0_R_WAKE`, `FIN_TIP_R` |
| cuerpo más largo | `UPSTREAM_L`/`DOWNSTREAM_L` son en longitudes de cuerpo: el dominio crece solo. `X_WAKE_2` es absoluto |
| aletas más largas en cuerda | `FIN_X_LEAD` (el bloque refinado tiene que empezar sobre el cilindro) |
| aletas más altas | una zona radial fuera de `FIN_TIP_R` (`--preset fintip` como modelo) |
| otro perfil de cono | reemplazar `r_nose`; mantener monotonía y `r_nose(L_NOSE) == R_BODY` |

## Otro cohete

Copiar `aconcaguaGeom.py` a `<nombre>Geom.py`, cambiar constantes y, si hace
falta, `r_nose`, y apuntar el `import aconcaguaGeom as G` de `meshParams.py`,
`buildHexBody.py`, `finPatch.py`, `meshFinish.py` y `build.py` al nuevo
módulo (es un `import`, no un parámetro, a propósito: la geometría de una
malla no es algo que se cambie por línea de comando). La topología sirve para
cualquier cuerpo de revolución con punta, cilindro, boattail o no, base plana
y cuatro aletas a 90°.

## Valores de referencia

`reference_values(fraction)`: `Aref = π D²/4 × fraction` (sección del cuerpo
por la fracción de 360° que hay en la malla), `lRef = D`. Convención de
cohetería (Barrowman, OpenRocket): **no** la sección de las aletas. Un `Aref`
anterior usaba `R = 0.2355` (la semi-envergadura) y daba `Cd` 9.7 veces bajo.
