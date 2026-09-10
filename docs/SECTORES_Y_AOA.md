# Sectores y ángulo de ataque

La estructura de bloques se construye siempre para **un cuadrante** (θ de 0° a 90°) con media aleta en cada uno de sus dos planos. Una malla de 180° o de 360° se arma **rotando copias** de ese cuadrante y **cosiéndolas** nodo a nodo en los planos donde se encuentran. Los parámetros, la calidad y los nombres de patch son idénticos en los tres casos; sólo cambian el número de celdas y
qué ángulos de viento se pueden simular.
## Ejes y convenciones

```
                 +y  (aleta 0)
                  │
                  │      cuadrante (sector quarter): y ≥ 0, z ≤ 0
   +z ────────────┼──────────── -z (aleta 1)     planos de simetría: z = 0 (A)  y  y = 0 (B)
   (aleta 3)      │
                  │
                 -y  (aleta 2)
```

- Eje del cuerpo +x, punta en x = 0, base en x = 2.955 m. El viento viene de −x.
- Las cuatro aletas están en +y, −z, −y, +z (a 0°, 90°, 180°, 270°).
- `alpha` rota la corriente libre de +x hacia +y (plano de cabeceo x–y, que
  contiene el par de aletas +y/−y): configuración **"+"**.
- `beta` la rota hacia +z. `alpha = beta` es un balanceo de 45°: configuración
  **"×"**. Para un ángulo total σ con balanceo φ respecto del plano +y:
  `tan α = tan σ · cos φ`, `sin β = sin σ · sin φ`.

## Los tres sectores

| `--sector` | Copias | Dominio | Planos `symm` | Permite | Celdas (malla fina) |
|---|---|---|---|---|---|
| `quarter` | 1 | y ≥ 0, z ≤ 0 | z = 0 y y = 0 | sólo flujo axial: `alpha = beta = 0` | 3.54 M |
| `half` | 2 | z ≤ 0 | z = 0 | `alpha ≠ 0`, `beta = 0` | 7.08 M |
| `full` | 4 | 360° | ninguno | cualquier `alpha`, `beta`; estela transitoria asimétrica | 14.1 M |

```bash
python build.py --preset medium --sector half     # 2.4 M celdas, ~2 min
python build.py --preset medium --sector full     # 4.9 M celdas, ~4 min
```

`Allrun` lee `sector` de `constant/meshInfo` y tira error si desde `flowConditions` se pide un ángulo que la malla no puede representar.

`Aref` en el `meshInfo` es la sección del cuerpo por ¼, ½ o 1, así que `Cd`, `Cl`, `Cm` salen en su valor verdadero para cualquier sector. Las fuerzas dimensionales (`forces1`) hay que multiplicarlas por 4, 2 o 1.

## Cómo se arma (sectorAssembly.py)

1. **Copias.** La copia k es el cuadrante rotado −90°·k alrededor de x (rotación propia: la orientación de los hexaedros se conserva).
2. **Correspondencia de nodos.** El plano B (y = 0) de la copia k cae sobre el plano A (z = 0) de la copia k+1. Dos nodos de esos planos son el mismo nodo cuando tienen la misma (x, r), y la tienen exactamente porque el cuadrante es simétrico especularmente respecto de 45° (los anchos y las gradaciones azimutales están construidos así; `AZ_BLOCK_GROWTH` y `AZ_FIN_H` se
   aplican desde ambos planos, y la relajación `AZ_RELAX` conserva la simetría porque interpola entre dos distribuciones que la tienen). El emparejamiento usa un KD-tree con tolerancia 1e-7 m y exige que sea uno a uno.
3. **Las aletas.** Un nodo que la deformación de la aleta movió fuera del plano **no se fusiona**: queda como dos nodos, uno a +t/2 y otro a −t/2, que son las dos caras de la aleta de espesor completo. Las caras del plano cuyos cuatro nodos se fusionaron pasan a ser interiores y se descartan de ambos lados; las caras con algún nodo sobre la aleta son pared `fins` de ambos
   lados. Misma regla que clasifica las caras `symm`/`fins` del cuadrante `finPatch.split_symm`.
4. **Auditoría.** La malla ensamblada se pasa por la misma auditoría que el cuadrante: si algún par de caras no se fusionó, aparece como cara de borde sin patch y `build.py` no escribe. Con `--no-audit` se salta este paso.

```
sector: 2 x quadrant -> 68,025 nodes, 61,536 cells; 1 interface(s), 2,955 nodes merged per interface, 2,752 faces made interior, 316 fin faces
```

## Qué tener en cuenta

- **Sólo 4 aletas a 90°.** Las aletas viven en los planos del cuadrante. Tres aletas necesitan un sector de 120°, que un núcleo butterfly cuadrado no puede teselar: sería otra topología de bloques, no otro ensamblado.
- **Borde de ataque / fuga.** El contorno del patch `fins` está cuantizado a una celda, así que la aleta se extiende hasta una celda de más en todo el perímetro: 401.8 cm² contra 376.0 exactos por semi-aleta en la fina, +6.9 %. La **normal** de la pared es correcta en todos lados, así que no es la escalera que corrompe el esfuerzo de corte. Converge de primer orden con `FIN_H_X` y `FIN_H_R`. Detalle y efecto sobre los coeficientes en [WORKFLOW.md](WORKFLOW.md#el-borde-en-escalera-del-patch-fins).
- **Malla completa sin plano de simetría y flujo axial.** Es la elección para una estela transitoria (`pimpleFoam`): un plano de simetría suprime los modos asimétricos del desprendimiento.
- **Costo.** El build es serie y escala con las celdas: la completa fina tarda ~5 min y pesa ~2 GB en ASCII. `gmshToFoam` sobre ese archivo necesita ~15 GB de RAM. Para barridos en ángulo, `coarse half` (1.06 M) es el punto de partida razonable y `medium half` (2.4 M) el de trabajo.
