# CLAUDE.md — Aconcagua CFD

Aerodinámica del cohete Aconcagua en OpenFOAM v2412. Malla con cfMesh, una
plantilla de caso con `rhoCentralFoam`, barrido de Mach en serie.

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

Todo el barrido es a α = 0. No pagues celdas por el vórtice de punta de aleta.

**Una malla por punto de Mach.** El espesor de la primera capa sale de la
extrusión contra la celda de superficie, no de un y+ objetivo, así que el y+
real se mueve con la velocidad. Medilo en cada corrida.

---

## La estructura

```
flightModel/mesh/         makeRocket.py, stl/, Allmesh, meshInfo, system/
flightModel/common/       Allrun, Allrefine
flightModel/case-central/ la plantilla, rhoCentralFoam para todo Mach
flightModel/newCase.sh    plantilla + malla + condiciones -> una corrida
flightModel/run.sh        sweep.txt en serie, sobre una malla
```

El único Python es `mesh/makeRocket.py`, que escribe el STL; sólo usa la
biblioteca estándar. La malla se define en `mesh/system/meshDict` y todo
ahí es un **nivel**: un entero que parte una celda al medio, local a un patch o
a una caja, sin acoplarse a nada más.

---

## Reglas

- **Medí antes de proponer.** `mesh/Allmesh <dir>` da celdas, volúmenes
  negativos, no-ortogonalidad y skewness en unos minutos. Un cambio sin número
  al lado no está terminado.
- **`checkMesh` con volúmenes negativos = malla inservible.** Los otros
  chequeos de `-allGeometry` que fallan hoy (determinante chico en las celdas
  de capa, face tets, celdas cóncavas) son de otra categoría y no bloquean.
- **Las aletas van con `nLayers 3`.** Más enreda la extrusión. Se midió con
  la aleta vieja de 6 mm y borde de ataque en filo; la actual tiene 11.7 mm
  y bordes con radio de 1 mm, y no se volvió a medir.
- **Los cantos de las aletas llevan radio.** Con filo, las capas de las dos
  caras chocan en el borde (skewness 25–33). Está en `docs/MESH.md` §3.3.
- Antes de cortar o cambiar algo, grepeá `docs/` además del código.
- Una tarea por sesión. Un commit por tarea.

## Trampas conocidas

- **La repo vive adentro de OneDrive y OneDrive revierte archivos.** Pasó: una
  restauración de sync dejó un archivo fuente en una versión anterior a dos
  commits, y el archivo igual parseaba, así que no se notaba leyéndolo. Después
  de cualquier interrupción, `git diff --stat` antes de confiar en el working
  tree; un archivo con cientos de líneas borradas que vos no borraste es esto.
  Se recupera con `git checkout HEAD -- <archivo>`. El mitigante real es
  commitear seguido.
- **No metas scripts de Python en heredocs de bash**: el `\n` dentro de un
  f-string se convierte en salto de línea real y rompe el archivo.
- **`runApplication` devuelve 0 aunque la aplicación muera** con FOAM FATAL
  ERROR, así que `set -e` no lo agarra. `mesh/Allmesh` revisa cada log por eso.
- Un `*/` dentro del banner de comentario de un dict de OpenFOAM cierra el
  comentario antes de tiempo y el parser lee basura.
- No edites un script de shell mientras corre: `sh` lo lee a medida que lo
  ejecuta.
