# case-pimple/ — plantilla a bajo Mach, `rhoPimpleFoam`

Transitoria y basada en presión: la acústica es implícita, así que el paso lo
fija el Courant convectivo (`maxCo 1` sobre |U|) y no |U| + c. A M 0.3 eso da
un paso ~10 veces mayor que `case-central` sobre la misma malla.

```bash
../newCase.sh ~/runs/p03 --mesh ~/meshes/coarse --Minf 0.3 --template case-pimple
```

- Sin rampa de Courant ni `tauCo`: arranca en `deltaT 1e-6` con paso
  adaptativo.
- `endTime 0.06` s: dos pasadas del flujo sobre el cohete a M 0.3.
- Condiciones de borde subsónicas de `case-central` (`freestreamPressure`,
  `freestreamVelocity`). No usar con `Minf > 1`.
- `fvOptions` (recorte de T) sí se aplica con este solver.

**Estado:** probada en la malla gruesa sin capas (`~/meshes/coarse`). Sobre
mallas con capas límite divergió en las esquinas del borde de fuga de las
aletas: raíz y puntera, donde la malla deja celdas astilla.
