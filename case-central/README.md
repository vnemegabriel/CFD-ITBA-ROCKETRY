# case-central/ — la plantilla, M 0.2 a 1.8

`rhoCentralFoam` con flujo central-upwind de Kurganov, `kOmegaSST` con wall
functions, aire como gas perfecto con viscosidad de Sutherland. Una sola
plantilla para todo el barrido; se copia con `../newCase.sh`.

**No está validada.** Los números no tienen respaldo hasta que se haga
[../docs/VALIDATION.md](../docs/VALIDATION.md), y los de M 0.2 y 0.4 son los
más sospechosos: ver [../docs/SOLVERS.md](../docs/SOLVERS.md), §2.2.

```bash
../newCase.sh ~/runs/m02 --Minf 0.2
../newCase.sh ~/runs/m18 --Minf 1.8
```

- **Lo que cambia con el Mach** es sólo la condición de borde de `p`:
  `0.orig/p` elige `freestreamPressure` o `fixedValue`/`zeroGradient`/
  `waveTransmissive` según `Minf > 1`. Todo lo demás es igual. Ver
  [../docs/SOLVERS.md](../docs/SOLVERS.md), §7.
- **Arranca con rampa de Courant**, calcada del OpenFOAM ToolChain: `maxCo` de
  0.01 a `final_Co` y upwind → van Leer, con los parámetros en `config.json`.
  Ver [../docs/SOLVERS.md](../docs/SOLVERS.md), §2.3.
- **Es cara.** Paso global adaptativo con `Euler`: el paso lo fija el límite
  acústico de la celda más chica. Con la malla actual eso son ~2e-10 s por
  unidad de Courant.
- **La malla se dimensiona para una velocidad**, y del punto más lento al más
  rápido el barrido cambia 9 veces. Leé `Uinf` y el y+ esperado del banner de
  `./Allrun`.
