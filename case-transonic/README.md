# case-transonic/ — plantilla compresible, 0.3 < M < 1.2

`rhoSimpleFoam` con `transonic yes` + `kOmegaSST` con wall functions, aire como
gas perfecto con viscosidad de Sutherland. Se copia a un directorio de corrida
con `../newCase.sh --regime trans` y se corre ahí.

**No está validada.** Arranca y corre; los números que dé no tienen respaldo
hasta que se haga [../docs/VALIDATION.md](../docs/VALIDATION.md).

La condición se da como **número de Mach**, no como velocidad:

```bash
../newCase.sh ~/runs/m09 malla.msh --regime trans --Minf 0.9 --alpha 2
```

`system/flowDerived` deriva de ahí la velocidad del sonido, el vector
velocidad, la densidad, la viscosidad por Sutherland y la turbulencia de
entrada. Cada fórmula, con su verificación numérica, está en
[../docs/SOLVERS.md](../docs/SOLVERS.md).

## Qué cambia respecto de la subsónica

| | subsónica | ésta |
|---|---|---|
| solver | `simpleFoam` | `rhoSimpleFoam`, `transonic yes` |
| `p` | cinemática, m²/s² | **Pa** |
| campos extra en `0.orig/` | — | `T`, `alphat` |
| propiedades | `transportProperties`, ν constante | `thermophysicalProperties`, Sutherland |
| inicialización | `potentialFoam` | corriente uniforme |
| guarda numérica | — | `constant/fvOptions`, clip de temperatura |

`constant/fvOptions` no es un modelo físico sino una protección: sin él la
primera iteración desde una corriente uniforme lleva la temperatura fuera de
rango y el solver muere con una excepción de punto flotante dentro de la
librería termodinámica. Después de una corrida hay que verificar que el
mínimo y el máximo de `T` estén estrictamente dentro de los límites; si la
solución se apoya en uno, el resultado lo está fijando ese archivo y no la
física.

## Qué mirar

El campo `MachNo` que escribe el `controlDict`. Donde cruza 1 hay una burbuja
supersónica; si esa burbuja llega al boattail, el problema ya no es
transónico suave y conviene pasar a `../case-supersonic/`.
