# CLAUDE.md — Aconcagua CFD

Malla estructurada all-hex del cohete Aconcagua (gmsh → OpenFOAM v2412) más tres
casos: `case-subsonic`, `case-transonic`, `case-supersonic`.

Escribí en **castellano rioplatense**, no en inglés traducido. La jerga CFD queda
en inglés (mesh, solver, y+, wall function, skewness).

En los archivos, seguí el idioma del que estás tocando: el código y los
comentarios de `mesh/` están en inglés, y `docs/` en castellano salvo
`docs/MESH_DESIGN.md`, que quedó en inglés.

---

## El objetivo

**Curva de Cd contra Mach: 0.2, 0.4, 0.8, 1.2, 1.8**, y algún punto más si hace
falta. Lo que tiene que resolver bien la malla, en orden:

1. **El choque**, para M ≥ 0.8. Choque de proa, choque del borde de ataque de la
   aleta, expansión sobre el boattail. Es la fuente de la drag rise.
2. **La estela de base no estacionaria** para M > 1. El drag de base es una
   fracción grande del Cd total de un cohete y no es estacionario: eso empuja a
   `--sector full` y a `presets/wake_unsteady.py` en los puntos supersónicos.
3. La fricción de pared, que depende de y+ — ver abajo, es la trampa grande.

**Lo que NO es prioridad:** el vórtice de punta de aleta a incidencia. Todo el
barrido es a α = 0. No pagues celdas por bajar el chaflán de punta más allá de
lo que ya da `fintip` (6 mm); ese gasto se justificaría sólo con incidencia.

### ⚠ `y1` depende de U: es una malla por punto de Mach

`meshParams.derived()` calcula `y1` con `U` y `NU` de `meshParams.py`, hoy
`U = 100 m/s`. Con una sola malla para todo el barrido, el y+ real queda:

| Mach | U [m/s] | y1 correcto [mm] | y+ real con la malla de U = 100 |
|---|---|---|---|
| 0.2 | 68 | 0.4293 | **23** |
| 0.4 | 136 | 0.2301 | 42 |
| 0.8 | 272 | 0.1233 | 79 |
| 1.2 | 408 | 0.0856 | 113 |
| 1.8 | 612 | 0.0594 | **163** |

Un factor 7 de dispersión, con el punto de M 0.2 metido en la capa buffer donde
las wall functions son peores. Ese es un sesgo en la fricción **que depende de
Mach**, o sea la forma exacta de error que arruina una curva de Cd.

Construí **una malla por punto**, con `--set U=<Uinf>`. El banner de
`case-*/Allrun` ya lo dice y `case-*/system/flowConditions` lo repite.

---

## Verificación

```bash
python mesh/check.py            # ~1 min, el bucle de trabajo
python mesh/check.py --full     # ~6 min, antes de commitear
python mesh/check.py --update   # regenera check.baseline.json
```

**No digas que algo está listo sin haber corrido `--full`.** Dos capas:

- Un **piso absoluto** que nunca se mueve: cero volúmenes negativos, ninguna cara
  compartida por 3+ celdas, caras de borde == caras en parches.
- Un **baseline** (`mesh/check.baseline.json`) con los números que el repo
  produce hoy. Un umbral fijo no habría detectado las dos regresiones que
  realmente pasaron: el parche de aleta saliendo 10–21 % más grande que el
  planform, y la no-ortogonalidad media pasando de 5.78° a 7.59°.

Si un número se movió y estaba previsto, corré `--update` y **commiteá el
baseline junto con el cambio que lo movió**. El diff del baseline es la
evidencia.

---

## Qué es contrato

**`docs/` y `mesh/presets/` son interfaz soportada.** No saques nada de ahí ni le
cambies el comportamiento sin preguntar. Si algo está documentado, tiene usuario,
aunque ningún script del repo lo invoque: `--no-audit` existe porque auditar 14 M
celdas pide ~10 GB, y `--out` porque escribir el `.msh` dentro de OneDrive
mientras sincroniza duplica el tiempo. Los dos están en `docs/TROUBLESHOOTING.md`
y ninguno aparece en un `Allrun`.

Los **defaults de `meshParams.py` sí se pueden mover**, con números medidos, el
doc actualizado y `check.py --full` corrido.

Las **constantes de `aconcaguaGeom.py` son geometría del vehículo, no parámetros
de malla.** No las toques sin preguntar, ni siquiera por micrones.

---

## Doctrina del repo

- `aconcaguaGeom.py` es la **geometría de registro**. Reemplaza al STL: las
  fórmulas son la cosa, no una aproximación de ella. Todo se deriva de cinco
  números; si una constante se puede derivar de otra, se deriva.
- **Nada se snapea ni se relaja.** El nodo se coloca sobre la superficie
  evaluando la fórmula. No hay iteración que pueda quedar a medio camino.
- **Validar fuerte y temprano.** `validate_params()` y `G.validate()` corren
  siempre y nombran el parámetro y el valor que necesitás. Un acople silencioso
  es un bug aunque cada fórmula sea correcta por separado.
- Los comentarios explican **por qué**, y suelen citar el número del fallo que
  motivó la decisión. Mantené ese registro cuando cambies algo.
- La malla es 100 % hexaédrica y conforme. Si algo pide romper eso, es un cambio
  de topología y se discute, no se hace.

---

## Decisiones ya tomadas — no re-litigar

- **Los bordes de la aleta son fronteras de bloque**, no algo que la grilla cruza
  (`finEdge.py`). Se llegó ahí porque `split_symm()` sólo puede clasificar una
  cara entera, y esa regla está forzada por el stitch de sectores: una cara del
  plano se fusiona sólo si sus cuatro nodos siguen en z = 0. Un test por
  centroide bajaría el sesgo de área pero rompe el ensamble.
- **Los 57.4° de no-ortogonalidad en el LE (leading edge) son la flecha.** No bajan mientras las
  líneas radiales sean círculos. Bajarlos requiere un collar en O alrededor del
  planform, o sea otro generador.
- **La punta de aleta no puede ser cuadrada.** La aleta es una deformación
  azimutal; una punta cuadrada necesita una cara a r constante con celdas de un
  solo lado, y una estructura de bloques conforme no la da. Lo único ajustable es
  la altura del chaflán = `max(FIN_TIP_SMEAR, celda radial en la punta)`.
- **`FIN_X_LEAD` está acoplado a la flecha**: `estiramiento = 1 + 0.2508 / L`.
  `validate_params()` lo chequea contra `FIN_LEAD_MAX_STRETCH`.
- **`gmsh Curve In Surface` no sirve acá.** Las entidades embebidas sólo andan con
  Delaunay y HXT; con transfinito son incompatibles.
- El refinamiento **de volumen** alrededor de la aleta no es trabajo del
  generador: es `case-*/Allrefine` (`topoSet` + `refineMesh` 2:1), que lee la
  geometría de `constant/meshInfo`.

**Dos defectos abiertos, diagnosticados y medidos**, en
`docs/MESH_DESIGN.md` § *Lo que queda abierto*: el tamaño de celda se invierte
al cruzar `FIN_TIP_R` sobre la nariz (salto 1.59 en el cap), y la distribución
azimutal de las aletas se aplica a todo el cuerpo (arco max/min = 14.4). Los
números ya están: no hace falta volver a medirlos.

---

## Cómo trabajar acá

- **Una tarea por sesión.** El costo es (tamaño del contexto) × (cantidad de
  requests); arrastrar tres tareas anteriores multiplica las dos cosas.
- Antes de cortar o cambiar algo, **grepeá `docs/` además del código**. Ahí se
  perdieron dos round-trips.
- Medí antes de proponer. Este repo tiene el bucle armado: `build.py --no-write`
  da celdas, no-ortogonalidad, skewness y área de aleta contra la analítica en
  17 s (smoke) o 60 s (fine). Un cambio sin número al lado no está terminado.
- Verificá con la malla más chica que detecte el problema. `coarse` alcanza para
  casi todo; `fine` es para la confirmación final.
- Un commit por tarea. No mezcles cuatro trabajos en un working tree.

### Trampas conocidas

- **No metas scripts de Python en heredocs de bash.** El `\n` dentro de un
  f-string se convierte en salto de línea real y rompe el archivo. Escribí el
  script con la herramienta de archivos y corrélo.
- `presets/smoke.py` (H_SCALE 6) es incompatible con un `FIN_H_R` chico: la zona
  entre la punta y `ZONE_R[0]` queda más angosta que la celda escalada.
  `check.py` corre smoke con `FIN_H_R=None` por eso.
- El `.msh` de `fine` pesa ~460 MB por cuadrante. Escribilo fuera de OneDrive con
  `--out`.
