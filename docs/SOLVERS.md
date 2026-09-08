
> **Estado.** HAY QUE VALIDAR TODO!!!!
# Solvers y modelos: los tres regímenes

Qué resuelve cada plantilla, con qué modelo, y de dónde sale cada cosa. El objetivo de este documento es que puedas **rehacer las cuentas a mano** y verificar que el caso hace lo que dice.

```bash
./newCase.sh ~/runs/x <malla.msh> --regime sub     # M < 0.3    simpleFoam
./newCase.sh ~/runs/x <malla.msh> --regime trans   # 0.3 - 1.2  rhoSimpleFoam
./newCase.sh ~/runs/x <malla.msh> --regime super   # M > 1.2    rhoCentralFoam
```

## 1. Panorama

|                    | `sub`                     | `trans`                            | `super`                              |
| ------------------ | ------------------------- | ---------------------------------- | ------------------------------------ |
| directorio         | `case-subsonic/`          | `case-transonic/`                  | `case-supersonic/`                   |
| solver             | `simpleFoam`              | `rhoSimpleFoam`                    | `rhoCentralFoam`                     |
| formulación        | incompresible             | compresible, basada en presión     | compresible, basada en densidad      |
| tiempo             | `steadyState`             | `steadyState`                      | `Euler` transitorio, paso adaptativo |
| ρ                  | constante, no se resuelve | variable, de la ecuación de estado | variable, **es la incógnita**        |
| `p`                | cinemática, m²/s²         | Pa                                 | Pa                                   |
| energía            | no se resuelve            | `e` sensible                       | `e` sensible, dentro de `rhoE`       |
| campos en `0/`     | `U p k omega nut`         | `+ T alphat`                       | `+ T alphat`                         |
| viscosidad         | ν constante               | Sutherland μ(T)                    | Sutherland μ(T)                      |
| turbulencia        | k-ω SST                   | k-ω SST                            | k-ω SST                              |
| se parametriza con | `Uinf`                    | `Minf`, `pInf`, `Tinf`             | `Minf`, `pInf`, `Tinf`               |
| inicialización     | `potentialFoam`           | corriente uniforme                 | corriente uniforme                   |
El límite de 0.3 es la convención usual: el error de despreciar la compresibilidad en la presión es del orden de M²/4, o sea 2.3 % a M 0.3. 
Anderson, *Fundamentals of Aerodynamics*, 6ª ed.  Sección 8.3, da el desarrollo `p₀/p = 1 + M²/4 + ...` del que sale ese número.

## 2. Ecuaciones y modelos, por régimen

## 2.0 Que són los factores de relajación 

Controlan y ajustan las actualizaciones de valor del siguiente paso en la simulación.

### 2.1 Subsónico: `simpleFoam`

Resuelve las RANS incompresibles estacionarias, con `p` dividida por la densidad, que es constante:

```
div(U) = 0
div(U⊗U) − div(νeff grad(U)) = −grad(p)          p ≡ presión / ρ
```

Acoplamiento presión-velocidad **SIMPLEC** (`consistent yes`), que es SIMPLE con la aproximación consistente de Van Doormaal y Raithby (1984), *Numer. Heat Transfer* 7, 147-163. Permite valores de relajación más altos.

**Viscosidad:** `transportModel Newtonian`, ν constante tomada de `flowConditions`. No hay temperatura, así que no hay ley de viscosidad.

### 2.2 Transónico: `rhoSimpleFoam`

RANS compresibles estacionarias, basadas en presión. En `fvSolution` se activa la condición transónica:

```
SIMPLE { transonic yes; }
```

Con eso la ecuación de presión gana un término convectivo, `div(phid,p)`, y cambia de carácter de elíptica a hiperbólica donde M > 1. Sin él, `rhoSimpleFoam` no pasa de M ≈ 0.7. La discretización de ese término está en `fvSchemes` como `div(phid,p) Gauss upwind`.

`pMinFactor` y `pMaxFactor` acotan la actualización de presión por iteración, para no tener un valor negativo que interrumpa la iteración.
### 2.3 Supersónico: `rhoCentralFoam`

**No resuelve una ecuación de presión.** Avanza las variables conservadas ρ, ρU y ρE con un flujo central-upwind, y después aplica implícitamente las correcciones difusivas. Por eso su `fvSolution` no tiene ni `SIMPLE` ni factores de relajación: un esquema explícito basado en densidad no tiene qué sub-relajar.

**Esquema de flujo:** `fluxScheme Kurganov`, el esquema central-upwind que parte el flujo de cara según las velocidades de onda locales `a⁺` y `a⁻`. Captura el choque sin resolver un problema de Riemann y sin viscosidad artificial ajustada a mano.

- Kurganov y Tadmor (2000), "New High-Resolution Central Schemes for Nonlinear Conservation Laws and Convection-Diffusion Equations", *J. Comput. Phys.* **160**, 241-282.
- Kurganov, Noelle y Petrova (2001), *SIAM J. Sci. Comput.* **23**(3), 707-740.
- La implementación en OpenFOAM y su validación: Greenshields, Weller, Gasparini y Reese (2010), "Implicit and explicit schemes for flows of aerodynamic and turbomachinery fluids", *Int. J. Numer. Meth. Fluids* **63**, 1-21.

**Reconstrucción:** las variables primitivas se llevan a las caras con limitador van Leer (`reconstruct(rho) vanLeer`, `reconstruct(U) vanLeerV`, `reconstruct(T) vanLeer`). El limitador baja el esquema a primer orden en la discontinuidad y lo mantiene de segundo orden en las zonas suaves, que es lo que hace que el choque quede monótono. van Leer (1979), *J. Comput. Phys.* **32**, 101-136.
## 3. Turbulencia: idéntica en los tres

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
| `alphat` | `compressible::alphatWallFunction`, `Prt 0.85` | solo en los casos compresibles: difusividad térmica turbulenta a partir de νt |

El objetivo es y+ ≈ 32 en el centro de la primera celda. **La malla se
dimensiona para una velocidad**: ver
[PARAMETERS.md](PARAMETROS.md#flujo-la-malla-se-dimensiona-para-una-velocidad).
`Allrun` estima el y+ que te va a quedar y avisa si se fue de banda.

## 4. Propiedades termofísicas (solo `trans` y `super`)

`constant/thermophysicalProperties` es idéntico en los dos casos compresibles:

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
| `hePsiThermo` | termodinámica basada en energía usando ψ = 1/(RT), la compresibilidad. Es lo que esperan tanto los solvers basados en presión como `rhoCentralFoam` |
| `pureMixture` | una sola especie, sin transporte de composición |
| `perfectGas` | `p = ρ R T` |
| `hConst` | Cp constante, `h = Cp(T − Tstd) + Hf` |
| `sutherland` | `μ(T) = As √T / (1 + Ts/T)` |
| `sensibleInternalEnergy` | la variable de energía resuelta es `e`, no `h`. `rhoCentralFoam` lo exige, y usarlo también en el transónico permite compartir este archivo |

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
| 0.8 | 325 K | < 0.5 % |
| 1.2 | 371 K | ~0.8 % |
| 1.8 | 475 K | ~1.8 % |
| 2.5 | 648 K | ~4 % |

Arriba de M 2.5 conviene cambiar `hConst` por `janaf`, que es lo que usan los
tutoriales de `rhoCentralFoam`. Los coeficientes JANAF para aire están en
`$FOAM_TUTORIALS/compressible/rhoCentralFoam/biconic25-55Run35/constant/thermophysicalProperties`.

`Prt = 0.85` es el número de Prandtl turbulento, en la condición de pared de
`alphat`. Es el valor estándar para aire; Kays (1994), *J. Heat Transfer*
**116**, 284-295, discute su variación real.

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

En el caso subsónico la parametrización es directa por `Uinf` y `nu`, sin
temperatura.

## 6. Esquemas de discretización

### Comunes a `sub` y `trans`

| Término | Esquema | Por qué |
|---|---|---|
| `ddt` | `steadyState` | se busca el estacionario |
| `div(phi,U)` | `bounded Gauss linearUpwind limited` | segundo orden sesgado a contracorriente. `bounded` resta `U·div(phi)`, que en un estacionario no convergido no es cero y desestabiliza. `limited` nombra el gradiente limitado por celda, que impide que la reconstrucción se pase del rango de los vecinos |
| `div(phi,k)`, `div(phi,omega)` | `bounded Gauss upwind` | primer orden a propósito: k y ω son estrictamente positivos y el modelo está calibrado con upwind. Un esquema de orden alto acá compra una precisión que el cierre no tiene |
| `laplacian` | `Gauss linear limited corrected 0.33` | corrección no ortogonal completa hasta un limitador de 0.33, que es lo indicado para una malla con algunos miles de caras arriba de 70° |
| `snGrad` | `limited corrected 0.33` | idem |

`nNonOrthogonalCorrectors 1` acompaña a esos limitadores. Esta malla tiene
7563 caras arriba de 70° de un total de 10.7 M, ver
[MESH_DESIGN.md](MESH_DESIGN.md#current-quality).

### Propios de `super`

Ver §2.3. `fluxScheme`, las `reconstruct(...)` y `localEuler` son las entradas
que hacen el trabajo; `divSchemes` casi no se usa porque el esquema es
basado en densidad.

## 7. Condiciones de borde

La diferencia grande entre regímenes está acá, y viene de la teoría de
características: **cuántas ondas entran y cuántas salen por cada borde**.

### `sub` y `trans`: bordes subsónicos

`inlet`, `outlet` y `box` usan el mismo par, porque en subsónico cada borde
puede ser entrada o salida según el ángulo de ataque y la posición:

| Campo | Condición |
|---|---|
| `U` | `freestreamVelocity` |
| `p` | `freestreamPressure` |
| `T`, `k`, `omega` | `inletOutlet` |

Las dos primeras son `inletOutlet` con el valor de corriente libre: fijan
donde el flujo entra y extrapolan donde sale. Una sola entrada sirve para
cualquier ángulo.

### `super`: bordes supersónicos

| Borde | Condición | Razón |
|---|---|---|
| `inlet` | todo `fixedValue` | a M > 1 **todas** las características entran, así que todo se impone y nada se extrapola |
| `outlet` | todo `zeroGradient` | todas salen: imponer algo sería sobre-especificar y reflejaría hacia aguas arriba |
| `box`, `p` | `waveTransmissive` | la velocidad **normal** a este borde es casi cero, así que no es un borde supersónico en esa dirección y las ondas lo cruzan en ambos sentidos. `waveTransmissive` advecta la presión hacia afuera a `u + c` y relaja hacia `fieldInf` en la longitud `lInf` |
| `box`, `U` y `T` | `inletOutlet` | |

`waveTransmissive` implementa la condición no reflectante de Poinsot y Lele
(1992), "Boundary Conditions for Direct Simulations of Compressible Viscous
Flows", *J. Comput. Phys.* **101**, 104-129.

**Esto importa acá.** A M 1.8 el ángulo de Mach es 33.7° y el cono desde la
punta alcanza el farfield en x = 17.6 m, adentro del dominio, que llega a
41.4 m. O sea que la onda **sale por el borde lateral**, no por el outlet. Una
`fixedValue` ahí la reflejaría hacia adentro.

Las paredes son `noSlip` para `U` y `zeroGradient` para `p` y `T` en todos los
regímenes. `T` adiabática: a M 1.8 la temperatura de recuperación ronda los
470 K, que el fuselaje no alcanza a seguir en 4 s de quemado, así que flujo de
calor nulo es el límite correcto y no una temperatura de pared fija.

## 8. Coeficientes de fuerza

En incompresible `p` es cinemática, así que `forceCoeffs` la multiplica por
una densidad de referencia: `rho rhoInf; rhoInf 1.225;`.

En compresible `p` ya está en Pa y la densidad es un campo, así que la
integración usa el campo: `rho rho;`. `rhoInf` sigue haciendo falta, para la
presión dinámica que normaliza los coeficientes, y en las plantillas
compresibles sale derivada de `pInf` y `Tinf`.

`Aref`, `lRef` y la lista de paredes vienen de `constant/meshInfo` en los tres
casos, así que siempre coinciden con la malla. Convención de referencia y
centro de presión: [WORKFLOW.md](WORKFLOW.md), sección 6.

## 9. De dónde salió cada configuración

Además de la literatura citada arriba, las dos plantillas compresibles siguen
tutoriales que vienen con OpenFOAM v2412, así que podés diferenciarlas contra
un caso que la distribución mantiene:

| Plantilla | Tutorial de referencia |
|---|---|
| `trans` | `$FOAM_TUTORIALS/compressible/rhoSimpleFoam/aerofoilNACA0012` — aerodinámica externa compresible estacionaria |
| `super` | `$FOAM_TUTORIALS/compressible/rhoCentralFoam/biconic25-55Run35` — cono biconico supersónico externo, con datos experimentales |

**`constant/fvOptions`** no sale de ningún tutorial: lo agregué porque sin él
las dos plantillas compresibles mueren en la primera iteración. Arrancando
desde una corriente uniforme, el transitorio inicial saca la temperatura de
rango, y con `perfectGas` eso da una velocidad del sonido imaginaria y una
excepción de punto flotante dentro de la librería termodinámica. El recorte a
150-1200 K lo evita. **No es un modelo físico**: después de converger hay que
verificar que el mínimo y el máximo de `T` estén estrictamente adentro. En la
prueba transónica quedaron en 277 y 312 K, con holgura.

Lo que cambié respecto de ellos, y por qué:

- Los limitadores de `laplacian` y `snGrad` subidos a `limited corrected 0.33`,
  por la no-ortogonalidad de esta malla.
- `transonic yes` en la plantilla transónica. El tutorial del perfil no lo usa
  porque corre más abajo en Mach.
- Turbulencia `kOmegaSST` en la supersónica; el tutorial del biconico es
  laminar, porque ese experimento lo es.
- Condiciones de pared rarificadas (`maxwellSlipU`, `smoluchowskiJumpT`) del
  biconico **no** se usan acá: ese experimento es de baja densidad y este
  cohete vuela a densidad de nivel del mar.

## 10. Lo que todavía falta

- **Validación.** Ninguna de las tres está contrastada contra datos. El plan
  está en [VALIDATION.md](VALIDATION.md).
- **Paso de tiempo local en el supersónico.** Ver §2.3: haría el barrido en
  Mach mucho más barato y hoy no funciona en este caso.
- **Malla adaptada al cono de Mach.** Las zonas radiales son cilindros; para
  resolver el choque hacen falta celdas alineadas con el cono.
- **Chorro de la tobera.** El disco de la base es pared en las tres
  plantillas. Modelar el escape pide partirlo en `nozzle` y `base` y una
  condición de presión y temperatura totales, y la presión de base con motor
  encendido es muy distinta de la apagada.
- **Malla por rango de Mach.** `y1` sale de una velocidad; ver §3.
