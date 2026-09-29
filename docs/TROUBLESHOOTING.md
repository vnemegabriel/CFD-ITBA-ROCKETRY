# Problemas conocidos

Qué significa cada error y qué hacer. Los mensajes están citados como salen.

## La malla (`mesh/Allmesh`)

**`checkMesh` reporta volúmenes negativos**
Casi siempre son las capas límite de las aletas. El borde de ataque es un filo
y la aleta tiene 6 mm de espesor: la extrusión se enreda contra sí misma.
Medido, mismo `meshDict` salvo `fins { nLayers }`:

| `nLayers` en `fins` | celdas | vol ≤ 0 | no-ortog máx | skew máx | arista mín |
|---|---|---|---|---|---|
| 10 | 1,211,476 | **58** | 174.8° | 229.6 | 2.15 µm |
| **3** | 862,462 | **0** | 78.6° | 14.6 | 12.2 µm |
| 0 | 762,714 | 0 | 72.9° | 13.6 | 31.4 µm |

Con 3 alcanza. Si tocás la geometría de la aleta y vuelven a aparecer, bajá a
2 antes de tocar cualquier otra cosa.

**Esta tabla está vieja.** Se midió cuando el `meshDict` nombraba los patches
`nosecone/body/boattail` y cfMesh no los encontraba, así que ojiva, cuerpo y
cola no tenían ni refinamiento ni capas. Con los nombres corregidos (2.3 M
celdas):

| `nLayers` en `fins` | celdas | vol ≤ 0 | no-ortog máx | skew máx |
|---|---|---|---|---|
| 10 | 2,299,408 | 0 | 75.6° | 14.9 |
| 3 | 2,299,408 | 0 | 75.6° | 14.9 |
| 3, otra corrida | 2,299,313 | **168** | 178.6° | 313.7 |
| 0 | 2,299,313 | **132** | 169.2° | 851.9 |

Que 3 y 10 den las mismas celdas es porque cfMesh unifica `nLayers` entre
patches vecinos: las aletas heredaban las 10 del cuerpo. Los volúmenes
negativos estaban todos en el borde de fuga, cerca de la puntera (x 2.88–2.93,
r 0.19–0.235). Se arregla con `allowDiscontinuity 1` en `fins`, que deja a las
aletas con sus 3 capas:

| `fins` | celdas | vol ≤ 0 | no-ortog máx | skew máx | caras con tets malos |
|---|---|---|---|---|---|
| `nLayers 3` (hereda 10) | 2,299,425 | **102** | 174.5° | 334.5 | 19,058 |
| `nLayers 3; allowDiscontinuity 1` | 2,075,194 | **0** | 87.5° | 14.6 | 3,492 |

Consecuencia: el primer espesor en las aletas es mayor que antes, y el y+
también. Ver la sección y+.

**`checkMesh` falla 6 chequeos pero sin volúmenes negativos**
Es lo normal en esta malla: skewness ~15, determinante chico en ~490 k celdas,
caras con peso de interpolación bajo. Son chequeos de `-allGeometry`. Lo que
hace inservible una malla son los volúmenes negativos y las caras compartidas
por tres celdas; eso está en cero.

**Un paso falla pero `Allmesh` sigue**
Ya no: `runApplication` devuelve 0 aunque la aplicación muera con
`FOAM FATAL ERROR`, así que `Allmesh` revisa cada log. Si ves un paso que
"pasó" y el siguiente se queja de un archivo que no existe, mirá el `log.` del
anterior.

**`FOAM FATAL IO ERROR: problem while reading header`**
Un `*/` dentro del banner de comentario de un dict cierra el comentario antes
de tiempo y OpenFOAM lee basura. Pasa al escribir rutas tipo `case-*/algo` en
el encabezado.

**La malla sale con el dominio equivocado**
Los seis números de `surfaceGenerateBoundingBox` son **metros medidos desde la
bounding box del modelo**, no múltiplos de ella. El utilitario imprime las dos
cajas: mirá `log.surfaceGenerateBoundingBox`.

**Quiero otra geometría**
Ponela en `mesh/stl/` y `mesh/Allmesh ~/meshes/otro mesh/stl/mi_cohete.stl`. El STL tiene que traer los solids `nosecone`,
`body`, `boattail` y `fins`; `Allmesh` los renombra a `cone`, `walls`, `tail`
y `fins`.

## El caso

**`newCase.sh` dice `no mesh in <dir>`**
Ese directorio no tiene `constant/polyMesh`. Las mallas no se versionan:
`mesh/Allmesh <dir>`.

**`Allmesh` dice `<dir> exists`**
Nunca pisa una malla. Elegí otro nombre o borrala.

**`newCase.sh` se niega a arrancar**
Nunca pisa un directorio existente. Elegí otro nombre o borralo.

**`!! <plantilla> has no flowConditions entry 'X'`**
Le pasaste una condición que la plantilla no usa. La velocidad es derivada:
se pasa `--Minf`, no `--Uinf`.

## Refinado local (`./Allrefine`)

**`!! constant/meshInfo carries no fin geometry (finTipR, finX0, ...)`**
El `constant/meshInfo` del caso es viejo. Copiá `mesh/meshInfo` encima.

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

## Casos compresibles

**El solver muere en la primera iteración con `Floating point exception` y una
traza que termina en `libfluidThermophysicalModels`**
La temperatura se fue de rango en el transitorio inicial desde la corriente
uniforme, y con `perfectGas` eso da una velocidad del sonido imaginaria. Lo
cubre `constant/fvOptions`, que recorta T entre 150 y 1200 K. Si el caso no lo
tiene, copialo de la plantilla. Si lo tiene y aun así muere, bajá los factores
de relajación de `fvSolution`.

**Después de converger, `min(T)` o `max(T)` coinciden con los límites de `fvOptions`**
Entonces el resultado lo está fijando ese archivo y no la física. Ampliá los
límites y volvé a correr; si la solución se sigue apoyando en ellos, hay algo
más mal.

**Con `ddtSchemes localEuler`, el Courant del log es absurdo, del orden de 10⁵**
Es esperado y no es un síntoma. Con paso local cada celda avanza con su propio
paso, y el Courant que imprime el log se calcula con el paso global, que ya no
gobierna nada. La plantilla viene con `Euler`.

**`rhoCentralFoam` muere con `sigFpe` en el primer paso**
El primer paso se da a `1.2 · deltaT` sin mirar `maxCo`. `deltaT` sale de
`tauCo` en `config.json`: medilo en esta malla ([SOLVERS.md §2.3](SOLVERS.md)).

## Lo que se ve en ParaView

**La aleta se ve gruesa o escalonada**
Subí el nivel de `fins` en `localRefinement` del `meshDict`. El borde de
ataque lo sostiene la extracción de aristas de feature
(`surfaceFeatureEdges -angle 30`), no el nivel: si el borde se ve redondeado,
el problema es el ángulo, no el refinamiento.

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

**y+ no lo fijás vos**
En cfMesh la capa límite se extruye: el primer espesor sale de `nLayers` y
`thicknessRatio` contra la celda de superficie local, no de un y+ objetivo.
Medilo en cada corrida con la función `yPlus` y ajustá `nLayers` o el nivel de
refinamiento del patch.

**y+ distinto en cada punto de Mach**
Esperable, y es la razón por la que el barrido usa una malla por punto: el
primer espesor es geométrico y el y+ escala con la velocidad. Si el barrido
tiene que comparar fricción entre puntos, ajustá `nLayers` por punto hasta que
el y+ quede en la misma banda.

**y+ alto en `fins` y razonable en el cuerpo**
Las aletas van con `nLayers 3` y el cuerpo con 10, así que el primer espesor
ahí es mayor. Subir las aletas por encima de 3 enreda la extrusión: subí en su
lugar el nivel de `fins` en `localRefinement`, que achica la celda de
superficie y con ella la primera capa.

## Entorno

**`newCase.sh`: `!! needs jinja2`**
La plantilla renderiza sus `*.j2` con `jinja2-cli`: `pipx install jinja2-cli`
(con `sudo apt install pipx` antes; Ubuntu no deja usar `pip` sobre el
Python del sistema, y `jinja2-cli` no está en apt).

**Scripts con `\r`: `bad interpreter` o `syntax error near unexpected token`**
Finales de línea CRLF. El repo lleva `.gitattributes` con `eol=lf`; si un
archivo se editó con una herramienta de Windows que los cambió: `sed -i 's/\r$//' archivo`
o `dos2unix`.

**Corridas lentas, `processor*` que no se borran, archivos bloqueados**
La corrida está en OneDrive o en `/mnt/c`. Mover a `~/runs` en WSL.
