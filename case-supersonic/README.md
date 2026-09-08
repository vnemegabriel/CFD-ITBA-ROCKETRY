# case-supersonic/ — plantilla compresible, M > 1.2

`rhoCentralFoam` con flujo central-upwind de Kurganov, `kOmegaSST` con wall
functions, aire como gas perfecto con viscosidad de Sutherland. Se copia con
`../newCase.sh --regime super`.

**No está validada.** Arranca, marcha y baja residuales; los números no tienen
respaldo hasta que se haga [../docs/VALIDATION.md](../docs/VALIDATION.md).

```bash
../newCase.sh ~/runs/m18 malla.msh --regime super --Minf 1.8
```

## En qué se diferencia de las otras dos

Es **basada en densidad**: no resuelve una ecuación de presión. Avanza ρ, ρU y
ρE con el flujo central-upwind y después aplica las correcciones difusivas.
Por eso `system/fvSolution` no tiene ni `SIMPLE` ni factores de relajación, y
`system/fvSchemes` tiene `fluxScheme` y `reconstruct(...)` en lugar de
`divSchemes` interesantes. Las referencias del esquema están en
[../docs/SOLVERS.md](../docs/SOLVERS.md).

> **No llega al estacionario todavía.** Sobre la malla de prueba `smoke half`
> arranca bien, mantiene el Courant en 0.30 y estaciona el `Cd` cerca de 0.78,
> pero se cae con una excepción de punto flotante a t = 8.6e-05 s, un 2 % del
> arranque. Probé que no es el paso de tiempo ni la condición de farfield.
> Lo primero a intentar es una malla mejor y dimensionada para el Mach:
> `python3 build.py --preset coarse --sector half --set U=612`. El detalle
> está en [../docs/SOLVERS.md](../docs/SOLVERS.md), §2.4.

**Es cara, y hay que saberlo antes de lanzarla.** Se entrega con `Euler` y
paso global adaptativo: una marcha transitoria hacia el estacionario. El paso
lo fija el límite acústico de la celda más chica, y llegar al estacionario
pide del orden de diez tiempos de paso, 48 ms a M 1.8. Son cientos de miles de
pasos.

El acelerador estándar es el paso de tiempo local (`ddtSchemes localEuler`),
que es lo que usa el tutorial del biconico. **Acá no funcionó**: no controla el
paso, el solver reporta Courant del orden de 10⁶ y diverge en menos de diez
iteraciones. Por eso no viene activado. Si lo hacés andar, es la mejora más
grande disponible; mientras tanto, para un barrido en Mach conviene aflojar el
espaciado de pared, que es lo que fija el paso.

## Dos cosas que esta plantilla NO arregla, y son de la malla

**El cono de Mach sale por el borde lateral.** A M 1.8 el ángulo de Mach es
33.7° y el cono desde la punta alcanza el farfield en x = 17.6 m, adentro del
dominio. La condición de borde es no reflectante (`waveTransmissive`), pero
las celdas ahí miden 1.6 m, así que la onda se difumina mucho antes de llegar.
Para arrastre alcanza; para la estructura del choque, no.

**La malla se dimensiona para una velocidad.** A M 1.8 son 612 m/s. Corriendo
eso en la malla de 100 m/s el y+ queda en 163 en vez de 32. Leé `Uinf` del
banner de `./Allrun` y reconstruí con `--set U=<esa velocidad>`.
