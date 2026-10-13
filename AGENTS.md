# AGENTS.md — Aconcagua CFD

Aerodinámica del cohete Aconcagua en OpenFOAM v2412. Malla con cfMesh, una
plantilla `rhoPimpleFoam` para M < 1 y una `rhoCentralFoam` para M > 1 que
arranca desde una corrida pimple al mismo Mach, barrido de Mach en serie.

Escribí en **castellano rioplatense**, no en inglés traducido. La jerga CFD
queda en inglés (mesh, solver, y+, wall function, skewness). Los `docs/` están
en castellano; los scripts, en inglés.

**Pocos comentarios en el código.** Para explicar están los docs. Un comentario
se justifica sólo cuando el código no lo puede decir — por ejemplo que los seis
números de `surfaceGenerateBoundingBox` son metros y no múltiplos de la
bounding box.

---

## El objetivo

**Curva de Cd contra Mach: 0.3, 0.6, 0.8, 1.2, 1.8.** Lo que tiene que
resolver bien la malla, en orden:

1. El choque, para M ≥ 0.8.
2. La estela de base, que es una fracción grande del Cd de un cohete.
3. La fricción de pared, que depende de y+.

El barrido va a α = 0, 8 y 16°, impuesto con el campo de velocidad
(`--alpha`): la malla es la misma. Los conos de estela del `meshDict` siguen
el eje x, y a 16° la estela sale de ellos a menos de 0.4 m de la base.
No pagues celdas por el vórtice de punta de aleta.

**Un solo juego de mallas para todo el barrido:** h1, h2, h3 para la
convergencia (sólo cambia `maxCellSize`, razón √2) y h0, más gruesa, para
inicializar las corridas de h1 con `newCase.sh --init`. Media malla, con plano
de simetría en z = 0 y aletas a 45°: `Allmesh` malla el cohete completo y lo
corta en z = 0 con `subsetMesh`. Las capas dejan la primera celda del cuerpo
igual en todas: t1 ≈ 0.10 mm medido (y ≈ 0.05 mm), y+ ~5 a M 0.3 y ~50 a
M 1.8; `nutUSpaldingWallFunction` cubre el buffer. Medí t1 en cada malla
nueva con `post/firstLayer.py`. No armes mallas por Mach. El y+ real medilo
en cada corrida.

---

## La estructura

```
flightModel/mesh/         stl/Aconcagua.stl, Allmesh, meshInfo, system/
flightModel/common/       Allrun, Allrefine
flightModel/case-central/ plantilla rhoCentralFoam, M > 1
flightModel/case-pimple/  plantilla rhoPimpleFoam, M < 1 y arranque de M > 1: sólo lo que difiere del central
flightModel/newCase.sh    plantilla + malla + condiciones -> una corrida
flightModel/runMesh.sh    meshes.txt: construye las mallas que faltan
flightModel/runCase.sh    un sweep*.txt en serie, sobre mallas ya hechas: sweepSupersonic.txt
                          en el socket 0 y sweepSubsonic.txt en el 1, a la vez
flightModel/post/         monitor.sh (gnuplot en vivo), sweepTable.sh (CSV),
                          firstLayer.py (t1 e y de la primera celda)
```

El Python es `mesh/makeRocket.py`, que escribe el STL, `mesh/snapPlane.py`,
que aplana el corte de la media malla, y `post/firstLayer.py`, que mide la
primera capa; sólo usan la biblioteca estándar. La
malla se define en `mesh/system/meshDict` y todo ahí es un **nivel**: un
entero que parte una celda al medio, local a un patch o a una caja, sin
acoplarse a nada más.

---

## Reglas

- **Medí antes de proponer.** `mesh/Allmesh <dir>` da celdas, volúmenes
  negativos, no-ortogonalidad y skewness en unos minutos. Un cambio sin número
  al lado no está terminado.
- **`checkMesh` con volúmenes negativos = malla inservible.** Los otros
  chequeos de `-allGeometry` que fallan hoy (skewness ~15, determinante chico
  en ~490 k celdas de capa) son de otra categoría y no bloquean.
- **Las aletas van con `nLayers 3`.** Más enreda la extrusión: el borde de
  ataque es un filo y la aleta tiene 6 mm de espesor. Está medido en
  `docs/TROUBLESHOOTING.md`.
- Antes de cortar o cambiar algo, grepeá `docs/` además del código.
- Una tarea por sesión. Un commit por tarea.

## Cuidados

- **Con `localEuler` el Courant del log no significa nada.** rhoCentralFoam lo
  calcula con el `deltaT 1` de `controlDict` (da ~1e7) y rhoPimpleFoam no lo
  escribe. Mirá `Flow time scale min/max`. Esa lectura equivocada dejó al
  central en `Euler` global hasta 2026-10-12. `docs/SOLVERS.md` §5.1.
- **No metas scripts de Python en heredocs de bash**: el `\n` dentro de un
  f-string se convierte en salto de línea real y rompe el archivo.
- **`runApplication` devuelve 0 aunque la aplicación muera** con FOAM FATAL
  ERROR, así que `set -e` no lo agarra. `mesh/Allmesh` revisa cada log por eso.
- Un `*/` dentro del banner de comentario de un dict de OpenFOAM cierra el
  comentario antes de tiempo y el parser lee basura.
- No edites un script de shell mientras corre: `sh` lo lee a medida que lo
  ejecuta.
