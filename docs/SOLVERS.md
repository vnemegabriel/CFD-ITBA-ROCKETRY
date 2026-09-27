
> **Estado.** HAY QUE VALIDAR TODO!!!!
# Solver y modelos

Un solo solver para todo el barrido, M 0.3 a 1.8: `rhoCentralFoam`, en la
plantilla `case-central/`. Este documento dice qué resuelve, con qué modelo y
de dónde sale cada cosa, para que puedas **rehacer las cuentas a mano** y
verificar que el caso hace lo que dice.

```bash
./newCase.sh ~/runs/m03 --Minf 0.3
./newCase.sh ~/runs/m18 --Minf 1.8
```

## 1. Panorama

| | |
|---|---|
| directorio | `case-central/` |
| solver | `rhoCentralFoam` |
| formulación | compresible, basada en densidad |
| tiempo | `Euler` transitorio, paso global adaptativo, con rampa de Courant (§2.3) |
| ρ | variable, **es la incógnita** |
| `p` | Pa |
| energía | `e` sensible, dentro de `rhoE` |
| campos en `0/` | `U p T k omega nut alphat` |
| viscosidad | Sutherland μ(T) |
| turbulencia | k-ω SST con wall functions |
| se parametriza con | `Minf`, `pInf`, `Tinf` |
| inicialización | corriente uniforme |

**Por qué uno solo.** Antes había tres plantillas (`simpleFoam`,
`rhoSimpleFoam`, `rhoCentralFoam`). Un Cd que sale de tres solvers distintos
cambia de discretización justo en los bordes del régimen, y el salto en la
curva no se puede separar de la física. Con un solo solver, la curva Cd–Mach
tiene la misma numérica en todos los puntos.

## 2. Ecuaciones y numérica

### 2.1 `rhoCentralFoam`

**No resuelve una ecuación de presión.** Avanza las variables conservadas ρ,
ρU y ρE con un flujo central-upwind y después aplica implícitamente las
correcciones difusivas. Por eso `fvSolution` no tiene ni `SIMPLE` ni factores
de relajación: un esquema explícito basado en densidad no tiene qué
sub-relajar. `rho`, `rhoU` y `rhoE` usan `diagonal`; sólo `U`, `e`, `k` y
`omega` tienen un solver lineal de verdad.

**Esquema de flujo:** `fluxScheme Kurganov`, el esquema central-upwind que
parte el flujo de cara según las velocidades de onda locales `a⁺` y `a⁻`.
Captura el choque sin resolver un problema de Riemann y sin viscosidad
artificial ajustada a mano.

- Kurganov y Tadmor (2000), "New High-Resolution Central Schemes for Nonlinear Conservation Laws and Convection-Diffusion Equations", *J. Comput. Phys.* **160**, 241-282.
- Kurganov, Noelle y Petrova (2001), *SIAM J. Sci. Comput.* **23**(3), 707-740.
- La implementación en OpenFOAM y su validación: Greenshields, Weller, Gasparini y Reese (2010), "Implicit and explicit schemes for flows of aerodynamic and turbomachinery fluids", *Int. J. Numer. Meth. Fluids* **63**, 1-21.

**Reconstrucción:** las variables primitivas se llevan a las caras con
limitador van Leer (`reconstruct(rho) vanLeer`, `reconstruct(U) vanLeerV`,
`reconstruct(T) vanLeer`). El limitador baja el esquema a primer orden en la
discontinuidad y lo mantiene de segundo orden en las zonas suaves, que es lo
que hace que el choque quede monótono. van Leer (1979), *J. Comput. Phys.*
**32**, 101-136.

### 2.2 El precio a bajo Mach

**Por qué el barrido arranca en M 0.3.** Es la hipótesis de compresibilidad:
el error de despreciarla en la presión es del orden de M²/4, 2.3 % a M 0.3
(Anderson, *Fundamentals of Aerodynamics*, 6ª ed., §8.3, del desarrollo
`p₀/p = 1 + M²/4 + ...`). Abajo de eso el Cd casi no cambia con el Mach, y
correrlo con un solver compresible sólo compra rigidez.

`rhoCentralFoam` está hecho para choques, y a M 0.3 y 0.6 se paga de dos
maneras. Hay que saberlo antes de mirar esos dos puntos de la curva.

**Rigidez.** El paso explícito lo fija la onda más rápida, `|U| + c`; lo que
interesa es el transporte a `|U|`. A M 0.3 cada paso avanza el flujo un
M/(1+M) ≈ 23 % de lo que avanzaría con un paso convectivo, y llegar al
estacionario cuesta del orden de 1/M veces más pasos que a M 1.

**Disipación.** Un esquema upwind o central-upwind sin precondicionar mete una
disipación que escala con `c` y no con `|U|`. Cuando M → 0, las fluctuaciones
de presión que genera son de orden M, mientras que las físicas son de orden
M². El resultado es un campo de presión suavizado de más. Guillard y Viozat
(1999), "On the behaviour of upwind schemes in the low Mach number limit",
*Computers & Fluids* **28**, 63-86.

En la práctica: los Cd a M 0.3 y 0.6 se validan contra datos o contra una
corrida de `rhoSimpleFoam` antes de confiar en ellos
([VALIDATION.md](VALIDATION.md)). Si la diferencia no es aceptable, lo
indicado es volver a un solver basado en presión sólo en ese tramo, con la
discontinuidad numérica que eso trae.

### 2.3 Arranque: rampa de Courant

Calcada del OpenFOAM ToolChain (`templates/rhoCentralFoam`): arrancar con
`maxCo` muy bajo y esquemas de primer orden, y subir de a pasos. La corriente
uniforme contra una pared no-slip es un transitorio violento, y es ahí donde
explota.

La hace el functionObject `timeActivatedFileUpdate`, en `system/fileUpdater`.
En cada paso mira el tiempo y, al cruzar un umbral, copia otro archivo encima
de `system/controlDict` o de `system/fvSchemes`; `runTimeModifiable` lo relee.

| archivo | qué es |
|---|---|
| `config.json` | los parámetros de la rampa; `newCase.sh` renderiza con él todos los `*.j2` |
| `system/controlDictBase.j2` | el `controlDict` de siempre, sin `maxCo` |
| `system/controlDict`, `fileUpdaterControlDict/controlDict_{1..7,final}.j2` | incluyen la base y agregan su `maxCo`: 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.5, recortados a `final_Co` |
| `system/fileUpdater.j2` | el calendario: qué archivo entra a qué tiempo |
| `system/fvSchemesUpwind`, `fvSchemesVanLeer` | reconstrucción de primer orden y la de van Leer |
| `system/FOs/FOrunTimeControl.j2` | corte por convergencia del promedio de `Cd`, sólo en `controlDict_final` |

El calendario es el del ToolChain, en **iteraciones**: rampa en 0, 20, 35, 50,
65, 80, 100; a `vanLeer` pasa a van Leer y repite la rampa desde 0.01 (en
+0, 20, 35, 150, 300, 800, 2000); a `minIter` entra `controlDict_final` con
`final_Co` y el `runTimeControl`.

**La diferencia con el ToolChain.** Allá corren `localEuler`, donde el tiempo
es la iteración. Acá es `Euler` en segundos reales, así que el template
convierte cada iteración en tiempo con `t = Σ nᵢ · maxCoᵢ · tauCo`. `tauCo` es
el paso que da Co = 1 en la peor celda; se lee del log como `deltaT / maxCo`
una vez que el Co se estabiliza. Si está mal, la rampa se corre en el tiempo,
pero el Co sigue acotado por `maxCo`. Con `tauCo` también se fija el `deltaT`
inicial, porque `rhoCentralFoam` da el primer paso a `1.2 · deltaT` sin mirar
el Co: con un `1e-8` fijo, esta malla arrancaba con Co = 61 y moría con
`sigFpe` en el paso 1.

| clave | valor | qué hace |
|---|---|---|
| `tauCo` | `2e-10` | s. Medido en la malla completa de 2.3 M a M 1.8. `c` domina `|U| + c`, así que a M 0.3 es del mismo orden |
| `vanLeer` | `2000` | iteración del cambio a van Leer; `-1` se queda en upwind |
| `minIter` | `3000` | iteración de `controlDict_final` |
| `final_Co` | `0.3` | techo de toda la rampa |
| `coeffs_fields`, `coeffs_variation`, `coeffs_range` | `Cd`, `0.001`, `300` | condición del `runTimeControl`. A α ≠ 0 se pueden sumar `Cl` y `CmPitch`; a α = 0 no, porque su promedio es ~0 y el criterio relativo no cierra nunca |

Medido: el Co queda clavado en 0.0100 y a las 20 iteraciones pasa a 0.0200.
El mismo log muestra el problema de fondo: **Co medio 5e-9 contra máximo
0.02**. Una sola región de celdas diminutas fija el paso de 2.3 M celdas;
con `tauCo = 2e-10` s, llegar al estacionario son ~10⁹ pasos. La rampa evita
la explosión, pero no hace barata la corrida: eso es de la malla, que hoy
tiene volúmenes negativos ([TROUBLESHOOTING.md](TROUBLESHOOTING.md)).

## 3. Turbulencia

```
simulationType  RAS;
RASModel        kOmegaSST;
turbulence      on;
```

**Modelo de cierre.** Hipótesis de Boussinesq: el tensor de Reynolds se modela como una viscosidad turbulenta por la tasa de deformación,

```
−⟨u'ᵢu'ⱼ⟩ = νt (∂Uᵢ/∂xⱼ + ∂Uⱼ/∂xᵢ) − (2/3) k δᵢⱼ
```

con `νt = a₁k / max(a₁ω, S F₂)`. Ese denominador es el **limitador de tensión cortante** que distingue al SST de un k-ω común, y es la razón de elegirlo para un cohete: acota la producción de νt en gradiente de presión adverso, que es exactamente lo que pasa en el boattail y en la unión aleta-cuerpo, donde un k-ε o un k-ω sin limitar predicen la separación tarde.

**Referencias.**

- Menter (1994), "Two-Equation Eddy-Viscosity Turbulence Models for Engineering Applications", *AIAA Journal* **32**(8), 1598-1605. El modelo original, con las funciones de mezcla F₁ y F₂.
- Menter, Kuntz y Langtry (2003), "Ten Years of Industrial Experience with the SST Turbulence Model", *Turbulence, Heat and Mass Transfer* **4**, 625-632. **Ésta es la versión que implementa OpenFOAM**, no la de 1994: cambian la constante de producción y el limitador. Si comparás contra literatura, fijate cuál de las dos usa.
- La implementación exacta, con todas las constantes, está en el código:
  `src/TurbulenceModels/turbulenceModels/RAS/kOmegaSST/`. Constantes por
  defecto: `alphaK1 0.85`, `alphaK2 1.0`, `alphaOmega1 0.5`, `alphaOmega2
  0.856`, `beta1 0.075`, `beta2 0.0828`, `betaStar 0.09`, `gamma1 5/9`,
  `gamma2 0.44`, `a1 0.31`, `b1 1.0`, `c1 10`, `F3 no`.

**Tratamiento de pared.** Funciones de pared en las cuatro paredes:

| Campo | Condición | Por qué |
|---|---|---|
| `nut` | `nutUSpaldingWallFunction` | la ley de Spalding es **continua en y+**, así que la primera celda puede caer en la subcapa, en la zona buffer o en la capa logarítmica sin un salto en τw. Spalding (1961), "A Single Formula for the Law of the Wall", *J. Appl. Mech.* **28**(3), 455-458 |
| `k` | `kLowReWallFunction` | también continua en y+, por consistencia con la anterior |
| `omega` | `omegaWallFunction` | mezcla las formas de subcapa y logarítmica de ω |
| `alphat` | `compressible::alphatWallFunction`, `Prt 0.85` | difusividad térmica turbulenta a partir de νt |

El objetivo es y+ ≈ 32 en el centro de la primera celda. **La malla se
dimensiona para una velocidad**, y entre M 0.3 y 1.8 la velocidad cambia 6
veces: ver [WORKFLOW.md](WORKFLOW.md#3-el-barrido). `Allrun` estima el y+ que
te va a quedar y avisa si se fue de banda.

## 4. Propiedades termofísicas

```
thermoType
{
    type            hePsiThermo;
    mixture         pureMixture;
    transport       sutherland;
    thermo          hConst;
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleInternalEnergy;
}
```

| Entrada | Qué significa |
|---|---|
| `hePsiThermo` | termodinámica basada en energía usando ψ = 1/(RT), la compresibilidad. `rhoCentralFoam` sólo acepta `psiThermo` |
| `pureMixture` | una sola especie, sin transporte de composición |
| `perfectGas` | `p = ρ R T` |
| `hConst` | Cp constante, `h = Cp(T − Tstd) + Hf` |
| `sutherland` | `μ(T) = As √T / (1 + Ts/T)` |
| `sensibleInternalEnergy` | la variable de energía resuelta es `e`, no `h`. `rhoCentralFoam` lo exige |

**Constantes y su verificación.**

| | valor | verificación |
|---|---|---|
| `molWeight` | 28.96 kg/kmol | R = 8314.46/28.96 = **287.05** J/(kg K) |
| `Cp` | 1005 J/(kg K) | γ = Cp/(Cp − R) = 1005/717.95 = **1.400** |
| `As` | 1.4792e-06 kg/(m s √K) | a 288.15 K: μ = 1.4792e-6·√288.15/(1+116/288.15) = **1.790e-05** Pa s, el valor estándar del aire |
| `Ts` | 116 K | |

La ley de Sutherland en su forma habitual es `μ = μref (T/Tref)^1.5
(Tref+S)/(T+S)` con μref = 1.716e-05 Pa s, Tref = 273.15 K y S = 110.4 K
(White, *Viscous Fluid Flow*, 3ª ed., ec. 1-36). La forma de dos parámetros de
OpenFOAM es algebraicamente la misma: `As √T/(1+Ts/T) = As T^1.5/(T+Ts)`.

**Límite de validez de `hConst`.** La temperatura de estancamiento es
`T₀ = T(1 + (γ−1)/2 · M²)`:

| M | T₀ | error de Cp |
|---|---|---|
| 0.3 | 293 K | despreciable |
| 0.8 | 325 K | < 0.5 % |
| 1.2 | 371 K | ~0.8 % |
| 1.8 | 475 K | ~1.8 % |
| 2.5 | 648 K | ~4 % |

En todo el barrido `hConst` alcanza. Arriba de M 2.5 conviene cambiarlo por
`janaf`, que es lo que usan los tutoriales de `rhoCentralFoam`; los
coeficientes para aire están en
`$FOAM_TUTORIALS/compressible/rhoCentralFoam/biconic25-55Run35/constant/thermophysicalProperties`.

`Prt = 0.85` es el número de Prandtl turbulento, en la condición de pared de
`alphat`. Es el valor estándar para aire; Kays (1994), *J. Heat Transfer*
**116**, 284-295, discute su variación real.

**`constant/fvOptions`** recorta T entre 150 y 1200 K. Arrancando desde una
corriente uniforme, el transitorio inicial puede sacar la temperatura de
rango, y con `perfectGas` eso da una velocidad del sonido imaginaria y una
excepción de punto flotante dentro de la librería termodinámica. **No es un
modelo físico**: después de converger, el mínimo y el máximo de `T` tienen
que quedar estrictamente adentro.

## 5. De qué se derivan las condiciones

`system/flowConditions` tiene **solo números**. `system/flowDerived` hace el
álgebra con `#eval`, y todo lo demás la incluye. Cada fórmula es verificable:

| Derivada | Fórmula | A M 0.8, ISA nivel del mar |
|---|---|---|
| `aInf` | `√(γ R T)` | 340.3 m/s |
| `Uinf` | `M · a` | 272.2 m/s |
| `rhoInf` | `p/(R T)` | 1.225 kg/m³ |
| `muInf` | Sutherland | 1.790e-05 Pa s |
| `nuInf` | `μ/ρ` | 1.461e-05 m²/s |
| `kInlet` | `1.5 (U·Ti)²` | intensidad isótropa |
| `nutInlet` | `nuRatio · ν` | |
| `omegaInlet` | `k/νt` | de la definición `νt = k/ω` |
| `alphatInlet` | `ρ νt / Prt` | |

## 6. Esquemas de discretización

Ver §2.1. `fluxScheme` y las `reconstruct(...)` son las entradas que hacen el
trabajo; `divSchemes` casi no se usa porque el esquema es basado en densidad.
`system/fvSchemes` no se edita: la rampa (§2.3) lo pisa con `fvSchemesUpwind`
al arrancar y con `fvSchemesVanLeer` después.

| Término | Esquema | Por qué |
|---|---|---|
| `ddt` | `Euler` | marcha transitoria con paso global; `localEuler` no funcionó en esta malla (§9) |
| `div(phi,k)`, `div(phi,omega)` | `Gauss upwind` | primer orden a propósito: k y ω son estrictamente positivos y el modelo está calibrado con upwind |
| `laplacian` | `Gauss linear limited corrected 0.33` | corrección no ortogonal hasta un limitador de 0.33, por la no-ortogonalidad de esta malla |
| `snGrad` | `limited corrected 0.33` | idem |

## 7. Condiciones de borde

Vienen de la teoría de características: **cuántas ondas entran y cuántas
salen por cada borde**. `U` y `T` se resuelven igual a cualquier Mach; `p` no,
y `0.orig/p` elige la rama con `#if ${{ $Minf > 1 }}`.

| Campo | `inlet`, `outlet`, `box` | Por qué sirve a cualquier Mach |
|---|---|---|
| `U` | `freestreamVelocity` | fija donde el flujo entra y extrapola donde sale |
| `T`, `k`, `omega` | `inletOutlet` | idem |

**`p`, M ≤ 1:** `freestreamPressure` en los tres bordes. Es lo contrario de
`U`: extrapola donde el flujo entra y fija donde sale. En subsónico entra una
onda acústica por la salida, así que la presión del outlet se impone.

**`p`, M > 1:**

| Borde | Condición | Razón |
|---|---|---|
| `inlet` | `fixedValue` | a M > 1 **todas** las características entran, así que todo se impone |
| `outlet` | `zeroGradient` | todas salen: imponer algo sería sobre-especificar |
| `box` | `waveTransmissive` | la velocidad **normal** a este borde es casi cero, así que las ondas lo cruzan en ambos sentidos. Advecta la presión hacia afuera a `u + c` y relaja hacia `fieldInf` en la longitud `lInf` |

`freestreamPressure` en v2412 no tiene modo supersónico: fijaría `p` en un
outlet supersónico. Por eso la rama.

`waveTransmissive` implementa la condición no reflectante de Poinsot y Lele
(1992), "Boundary Conditions for Direct Simulations of Compressible Viscous
Flows", *J. Comput. Phys.* **101**, 104-129.

**Esto importa acá.** A M 1.8 el ángulo de Mach es 33.7° y el cono desde la
punta alcanza el farfield en x = 17.6 m, adentro del dominio, que llega a
41.4 m. O sea que la onda **sale por el borde lateral**, no por el outlet. Una
`fixedValue` ahí la reflejaría hacia adentro.

Las paredes son `noSlip` para `U` y `zeroGradient` para `p` y `T`. `T`
adiabática: a M 1.8 la temperatura de recuperación ronda los 470 K, que el
fuselaje no alcanza a seguir en 4 s de quemado, así que flujo de calor nulo es
el límite correcto y no una temperatura de pared fija.

## 8. Coeficientes de fuerza

`p` está en Pa y la densidad es un campo, así que la integración usa el
campo: `rho rho;`. `rhoInf` sigue haciendo falta, para la presión dinámica
que normaliza los coeficientes, y sale derivada de `pInf` y `Tinf`.

`Aref`, `lRef` y la lista de paredes vienen de `constant/meshInfo`, así que
siempre coinciden con la malla. Convención de referencia y centro de presión:
[WORKFLOW.md](WORKFLOW.md), sección 6.

## 9. De dónde salió la configuración

La plantilla sigue `$FOAM_TUTORIALS/compressible/rhoCentralFoam/biconic25-55Run35`
(cono biconico supersónico externo, con datos experimentales) y la rampa del
OpenFOAM ToolChain. Lo que cambia respecto del tutorial, y por qué:

- Limitadores de `laplacian` y `snGrad` en `limited corrected 0.33`, por la
  no-ortogonalidad de esta malla.
- Turbulencia `kOmegaSST`; el tutorial del biconico es laminar, porque ese
  experimento lo es.
- Las condiciones de pared rarificadas del biconico (`maxwellSlipU`,
  `smoluchowskiJumpT`) **no** se usan: ese experimento es de baja densidad y
  este cohete vuela a densidad de nivel del mar.
- `ddt Euler` en vez de `localEuler`. El paso local es el acelerador estándar,
  pero en esta malla no controló el paso: el solver reportó Courant del orden
  de 10⁶ y divergió en menos de diez iteraciones. Esa prueba fue sin rampa y
  sobre la malla con volúmenes negativos, así que vale repetirla.
- `fvOptions` con el recorte de T (§4).

## 10. Lo que todavía falta

- **Validación.** Nada está contrastado contra datos, y los puntos de bajo
  Mach son los más sospechosos (§2.2). El plan está en
  [VALIDATION.md](VALIDATION.md).
- **Paso de tiempo local.** Ver §9: haría el barrido mucho más barato,
  sobre todo a bajo Mach.
- **Malla sin volúmenes negativos.** Hoy fija un `tauCo` de 2e-10 s.
- **Malla adaptada al cono de Mach.** Las zonas radiales son cilindros; para
  resolver el choque hacen falta celdas alineadas con el cono.
- **Chorro de la tobera.** El disco de la base es pared. Modelar el escape
  pide partirlo en `nozzle` y `base` y una condición de presión y temperatura
  totales, y la presión de base con motor encendido es muy distinta de la
  apagada.
- **Malla por rango de Mach.** `y1` sale de una velocidad; ver §3.
