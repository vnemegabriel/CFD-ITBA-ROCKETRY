# case/ — plantilla OpenFOAM v2412

`simpleFoam` + `kOmegaSST` con wall functions, aire incompresible. Es una
**plantilla**: se copia a un directorio de corrida (`../newCase.sh`) y se
corre ahí. Flujo completo en [../docs/WORKFLOW.md](../docs/WORKFLOW.md).

```
case/
├── Allmesh   <msh>      gmshToFoam + meshInfo + fixPatchTypes + checkMesh + renumberMesh
├── Allrefine [tip|fins] OPCIONAL: refinado local 2:1 alrededor de las aletas
├── Allrun    [N]        chequeo sector/ángulos + decomposePar + potentialFoam + simpleFoam + reconstructPar
├── Allclean             vuelve a la plantilla
├── fixPatchTypes.py     pone wall / symmetry según constant/meshInfo.  NO opcional
├── 0.orig/              U p k omega nut  — todos #include system/flowDerived
├── constant/
│   ├── meshInfo         escrito por mesh/build.py, copiado por Allmesh: sector, Aref, patches
│   ├── transportProperties   nu de flowConditions
│   └── turbulenceProperties  kOmegaSST
└── system/
    ├── flowConditions   Uinf alpha beta nu rhoInf Ti nuRatio  — LO QUE SE EDITA
    ├── flowDerived      (Ux Uy Uz), drag/lift dirs, k omega nut  — #eval, no se edita
    ├── controlDict      simpleFoam, forceCoeffs/forces/yPlus/Cp/wallShearStress/residuals
    ├── fvSchemes        steady, linearUpwind, limited corrected 0.33
    ├── fvSolution       GAMG p, SIMPLEC, residualControl
    ├── topoSetDict      región a refinar; lee la geometría de las aletas del meshInfo
    ├── refineMeshDict   cómo refinarla (2:1, tres direcciones)
    └── decomposeParDict scotch, 8
```

## Refinado local (opcional)

`./Allrefine` corre **entre `Allmesh` y `Allrun`**, sobre una malla sin campos.
Refina 2:1 las celdas de un cilindro alrededor de las aletas, cuya posición
sale de `constant/meshInfo`, así que sigue a la geometría. Detalle y costos en
[../docs/WORKFLOW.md](../docs/WORKFLOW.md). Deja de ser 100 % hexaédrica: menos
del 1 % de las celdas pasan a poliedros en la transición, y la calidad medida
(no-ortogonalidad y skewness máximos) no cambia.

## Patches

| Patch | Tipo | Qué es | Existe en |
|---|---|---|---|
| `inlet` | patch | plano x = −17.7 m, `fixedValue` U | todos |
| `outlet` | patch | plano x = +41.4 m, `fixedValue` p | todos |
| `box` | patch | cilindro farfield r = 11.8 m, `freestreamVelocity`/`freestreamPressure` (inflow o outflow según el ángulo) | todos |
| `symm` | symmetry | planos z = 0 y/o y = 0 | quarter, half |
| `cone` | wall | ogiva | todos |
| `walls` | wall | cilindro | todos |
| `tail` | wall | boattail + disco de la base | todos |
| `fins` | wall | aletas | todos con `--fins` |

Las entradas `symm` en los campos se ignoran cuando el patch no existe; el
grupo `"(cone|walls|tail|fins)"` cubre las paredes con o sin aletas.

## Función objects

- `forceCoeffs1`: `Cd Cl Cm` con `Aref`/`lRef`/`patches` de `meshInfo` y
  direcciones de `flowDerived`. Momentos respecto de la punta.
- `forces1`: fuerzas en N (multiplicar por 4/2/1 según sector).
- `yPlus`: por patch, en el log y como campo.
- `Cp` (`pressure`, `staticCoeff`), `wallShearStress`: campos en cada escritura.
- `residuals` (`solverInfo`): `postProcessing/residuals/0/solverInfo.dat`.

## Cambiar de solver

La malla es agnóstica. Transitorio: `pimpleFoam`, `ddtSchemes backward`,
`deltaT ~ 5e-5` (Co ≈ 1 en la estela cercana fina), quitar `bounded` de los
`divSchemes`. Compresible: `rhoSimpleFoam`/`rhoPimpleFoam`, `p` en Pa,
agregar `T` y `thermophysicalProperties`; `flowConditions` sigue valiendo.
