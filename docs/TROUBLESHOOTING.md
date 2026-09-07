# Problemas conocidos

Qué significa cada error y qué hacer. Los mensajes están citados como salen.

## Construcción de la malla

**`refinement zones are inconsistent:` seguido de una lista**
`validate_params()` encontró parámetros incompatibles antes de construir.
Cada línea dice qué regla se rompió y con qué valores; la tabla de reglas
está en [PARAMETERS.md](PARAMETERS.md#validación). El caso típico es agregar
una zona a `ZONE_R` y olvidarse de `WAKE_ZONE_K`:
`ZONE0_R_WAKE 0.7 >= ZONE_R[1] 0.35: zone 0 would swallow zone 1 at the outlet`.

**`KeyError: presets/x.py: not a mesh parameter: ['ZONE_RR']`**
Un nombre mal escrito en un preset o en `--set`. Antes era un no-op
silencioso (la malla se construía con los defaults). La lista de nombres
válidos está en el mismo mensaje y en `meshParams.OVERRIDABLE`.

**`cell size (...) is not smaller than the segment it has to fill`**
Un `h_start`/`h_end` de `SEGMENTS` o un `ZONE_H` mayor que el tramo que tiene
que llenar, casi siempre por un `H_SCALE` grande sobre un tramo corto (`tail`
mide 125 mm). Bajar el tamaño o el `H_SCALE`.

**`interface nodes have no partner`** o **`interface planes carry N and M nodes`**
El ensamblado de sectores necesita que el cuadrante sea simétrico respecto de
45°. Aparece sólo tras editar `az_angles()`/`az_coefs()` en `meshParams.py`
o `_counts()` en `buildHexBody.py` de forma asimétrica. Si `scipy` no está
instalado el emparejamiento cae a redondeo exacto y puede fallar por ruido de
punto flotante: `pip install scipy`.

**`!! boundary faces and patch faces disagree -- refusing to write`**
Una cara de borde sin patch (iría a `defaultFaces` en OpenFOAM). En el
cuadrante no debería pasar nunca; en un ensamblado indica que un par de caras
de interfaz no se fusionó. Reportar con el `.params.py`.

**`ImportError: libGLU.so.1`** al importar gmsh en Linux/WSL:
`sudo apt install libglu1-mesa`.

**`MemoryError` o el proceso muere en la auditoría de la malla completa**
La auditoría de 14 M celdas necesita ~10 GB. `--no-audit` (el cuadrante se
audita igual y las copias tienen su misma calidad).

**El build tarda mucho más que lo tabulado**
`build.py` corre en serie. Con `mesh/output` en OneDrive, escribir 500 MB
mientras OneDrive los sube puede duplicar el tiempo: `--out` a un directorio
no sincronizado o excluir `output/` en OneDrive.

## Conversión y caso

**`gmshToFoam`: `Can only read ascii msh files`** o falla en `$MeshFormat`
El `.msh` no lo escribió `build.py` (que emite v2.2 ASCII). `gmsh` por
defecto escribe v4.1, que `gmshToFoam` no lee.

**`fixPatchTypes.py`: `!! UNRECOGNISED PATCHES: ['defaultFaces']`**
Caras sin grupo físico en el `.msh`. No seguir: reconstruir la malla.

**`fixPatchTypes.py`: `!! EXPECTED BUT ABSENT: ['symm']`** o `UNRECOGNISED PATCHES: ['symm']`
El `constant/meshInfo` no corresponde a esa malla (por ejemplo, un `meshInfo`
de un cuarto con un `.msh` completo). `Allmesh` copia el sidecar que está
**junto al `.msh`**; si se movió el `.msh` sin su `.meshInfo`, copiarlo a mano.

**`Allrun`: `!! quarter mesh: two symmetry planes, alpha and beta must be 0`**
Lo que dice. `--sector half` para `alpha`, `--sector full` para `beta`.
Ver [SECTORS_AND_AOA.md](SECTORS_AND_AOA.md).

**`forceCoeffs`: `Unknown patch name fins`** o `Cannot find patchField entry for symm`
`constant/meshInfo` y la malla no coinciden (ver arriba), o alguien editó
`wallPatches` a mano. En los campos de `0.orig/` la entrada `symm` puede
sobrar sin problema (OpenFOAM ignora entradas de patches que no existen),
pero **faltar** un patch en un campo sí es error.

**Cambié `alpha` y la corrida sigue con `U = (100 0 0)`**
Se usó `foamDictionary -set` sobre un archivo con `#eval` o `#include`
(`controlDict`, un campo, o `flowDerived`). `foamDictionary` reescribe el
archivo con todo evaluado y expandido, congelando los valores. Sólo
`system/flowConditions` (números planos) se edita con `foamDictionary`; el
resto con un editor o `sed`. Restaurar el archivo desde `case/` (o `git`).

**`checkMesh`: `Failed 3 mesh checks`**
Con `-allGeometry` es lo esperado en esta malla: aspect ratio > 1000 en
~3 000 celdas del farfield detrás de la base, determinante < 0.001 en las
celdas estiradas, y ~500 caras con peso de interpolación < 0.05. Ninguno
afecta a `simpleFoam`. Lo que sí debe cumplirse: 100 % hexaedros,
`Boundary openness OK`, no-ortogonalidad máx < 90°, skewness < 4. Detalle en
[WORKFLOW.md](WORKFLOW.md#qué-es-normal-en-checkmesh).

**`potentialFoam` o `simpleFoam` mueren en la primera iteración con `Phi` / `p`**
Casi siempre un patch de pared que quedó `type patch` (se saltó
`fixPatchTypes.py`) o el `symm` con `type patch`. `grep -A3 -E "^\s+(symm|cone|walls|tail|fins)" constant/polyMesh/boundary`.

**Diverge o la continuidad crece**
Bajar `relaxationFactors` (0.9 → 0.7 para `U`, 0.5 para `p` si se saca
`consistent`), o correr 200 iteraciones con `div(phi,U) bounded Gauss upwind`
antes de volver a `linearUpwind`. Con `alpha` grande (> 10°) el `potentialFoam`
inicial con `-initialiseUBCs` es importante; verificar que corrió (`log.potentialFoam`).

## Refinado local (`./Allrefine`)

**`!! constant/meshInfo carries no fin geometry (finTipR, finX0, ...)`**
El `.meshInfo` es anterior al soporte de refinado. Reconstruí la malla con
`mesh/build.py` para que el sidecar se escriba de nuevo; no hace falta cambiar
ningún parámetro.

**`!! there is a 0/ directory: refine BEFORE running`**
`refineMesh` mapea los campos que encuentra, así que corriendo después de
`Allrun` refinaría un resultado a medias en vez de la malla limpia. `./Allclean`
y empezar de nuevo, o refinar antes de la primera corrida.

**`!! constant/meshInfo says this mesh has no fins`**
La malla se construyó con `--no-fins`. No hay nada que refinar ahí.

**`checkMesh` ahora reporta poliedros**
Es lo esperado: la transición 2:1 no se puede representar con hexaedros. Son
menos del 1 % de las celdas y viven en la **superficie** del set, por eso un
anillo delgado genera casi tantos como una región mucho más grande. La
no-ortogonalidad y el skewness máximos no cambian.

**Después de refinar, el y+ bajo las aletas es la mitad que en el resto**
Usaste `./Allrefine fins`, cuya región llega hasta la pared del cuerpo. Para
una descomposición de arrastre por componente eso significa comparar patches a
dos y+ distintos. `./Allrefine tip` no toca la capa límite.

## Lo que se ve en ParaView

**La aleta se ve como un bulto borroso / escalonada (preset `coarse` o `--scale` > 1)**
Antes del 2026-09-06 los presets escalaban también la aleta: a `H_SCALE 3`
quedaba con 15 mm sobre la cuerda, una celda de espesor y 40 mm radiales en
la punta. Ahora `coarse` y `medium` llevan `H_SCALE_FIN = False` y
`FIN_H_R = 0.012`, que mantienen la aleta a la resolución de `fine` mientras
el resto se engrosa. Con `--scale` a mano, agregar `--set H_SCALE_FIN=False
--set FIN_H_R=0.012`. Ver [PARAMETERS.md](PARAMETERS.md#grosor-global).

**Celdas en "X" (moño) en un Clip con plano**
Aparecen cuando el plano del `Clip` pasa exactamente por nodos de la malla,
por ejemplo si el origen se eligió con `P` sobre un punto de la malla, o si
es un plano `y = const` tangente a una capa cilíndrica de nodos `r = const`.
ParaView tetraedriza los hexaedros que corta y esos slivers se dibujan como
moños. No son celdas de la malla: `checkMesh -allGeometry` da
`Face flatness min 0.99`, `Concave cell check OK`, `Face pyramids OK`.
Para mirar la malla usar `Clip` con **Crinkle clip** activado (muestra las
celdas enteras), o `Slice` con **Triangulate the slice** desactivado, o mover
el origen del plano fuera de una capa de nodos.

## y+

**y+ muy alto en `cone`/`walls`/`tail` con un preset**
`smoke` escala la pared también (`H_SCALE_WALL = True`): y+ ≈ 200. Es un
preset de prueba de pipeline, no de resultados. `coarse` y `medium` mantienen
`y1`.

**y+ muy alto en `fins` y razonable en el cuerpo**
En las aletas el primer espaciamiento normal lo fija `AZ_FIN_H`, no
`YPLUS_TARGET`. Con 0.5 mm da y+ ≈ 50 en la malla fina. Con `--scale` a mano
y `H_SCALE_FIN = True` escala con todo (1 mm y ~100 a `H_SCALE 2`); los
presets `coarse`/`medium` lo mantienen. Bajar `AZ_FIN_H` o subir `N_AZ_CELLS`.

## Entorno

**Windows: `python` vs `python3`**
En Windows el Python nativo es `python`; en WSL es `python3`. Los scripts de
shell (`Allmesh`, `newCase.sh`) llaman `python3` porque corren en WSL.

**Scripts con `\r`: `bad interpreter` o `syntax error near unexpected token`**
Finales de línea CRLF. El repo lleva `.gitattributes` con `eol=lf`; si un
archivo se editó con una herramienta de Windows que los cambió: `sed -i 's/\r$//' archivo`
o `dos2unix`.

**Corridas lentas, `processor*` que no se borran, archivos bloqueados**
La corrida está en OneDrive o en `/mnt/c`. Mover a `~/runs` en WSL.
