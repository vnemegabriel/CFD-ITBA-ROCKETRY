# Guía de Fabricación — Aletas 6 mm | Cohete Aconcagua

> **Para el equipo de manufactura.**  
> Este documento explica qué calculan los scripts de Python, cómo correrlos y cómo leer los resultados para ejecutar el proceso de laminado manual al vacío.

---

## Contenido

1. [Para qué sirven estos scripts](#1-para-qué-sirven-estos-scripts)
2. [Materiales y secuencia de laminado](#2-materiales-y-secuencia-de-laminado)
3. [Cómo correr los scripts](#3-cómo-correr-los-scripts)
4. [Cómo leer los resultados](#4-cómo-leer-los-resultados)
5. [Flujo completo de trabajo](#5-flujo-completo-de-trabajo)
6. [Parámetros que pueden ajustar](#6-parámetros-que-pueden-ajustar)
7. [Glosario](#7-glosario)

---

## 1. Para qué sirven estos scripts

Hay dos scripts principales. Deben correrse **en orden**:

### `inplaneG_v5.py` — Propiedades del laminado

Calcula cuán rígida es la aleta en función de cómo se apilan las capas de fibra de carbono. El resultado más importante para ustedes es la **rigidez torsional D66** y el **espesor total**.

> **En palabras simples:** dado el apilado de telas que eligieron, este script dice qué tan resistente a la torsión es la aleta.

### `flutterEstimate_v3.py` — Límite de flameo aeroelástico

Usa la rigidez calculada por el primer script para verificar que la aleta **no va a flamear** (vibrar de forma destructiva) durante el vuelo. Entrega un número llamado **Factor de Seguridad al Flameo (FSF)**.

> **En palabras simples:** este script dice si la aleta aguanta las cargas aerodinámicas del vuelo. FSF ≥ 1.50 significa que pasa. FSF < 1.50 significa que hay que rediseñar.

---

## 2. Materiales y secuencia de laminado

### Telas utilizadas

| Código interno | Descripción | Peso areal | Orientación |
|---|---|---|---|
| **CARBONODB300** | NCF biaxial carbono T700 | 300 g/m² | ±45° |
| **CARBONOGA90R** | Tejido plano carbono T700 | 302 g/m² | 0°/90° |

- Resina: **Epoxy LY1564**
- Fibra: **T700 12K**

### Secuencia de apilado — diseño AR1 (para aletas de 6 mm)

El laminado es **simétrico** respecto al plano medio. Se apila primero la mitad superior y luego se espeja para la mitad inferior.

**Media pila (mitad superior), de afuera hacia el núcleo:**

```
Capa  1 — DB300 ±45 NCF     (+45°)
Capa  2 — DB300 ±45 NCF     (−45°)
Capa  3 — GA90R tejido       (0°/90°)
Capa  4 — GA90R tejido       (0°/90°)
Capa  5 — DB300 ±45 NCF     (+45°)
Capa  6 — DB300 ±45 NCF     (−45°)
Capa  7 — GA90R tejido       (0°/90°)
Capa  8 — DB300 ±45 NCF     (+45°)
Capa  9 — DB300 ±45 NCF     (−45°)
Capa 10 — GA90R tejido       (0°/90°)
Capa 11 — GA90R tejido       (0°/90°)
──────── PLANO MEDIO ────────
[Espejo simétrico de las 11 capas anteriores]
```

**Total: 22 capas.** El laminado resultante a Vf = 0.50 tiene un espesor nominal de ~5.4 mm; escalando proporcionalmente se llega a **6.0 mm objetivo**.

### Fracción volumétrica de fibra (Vf)

El parámetro `Vf` representa cuánta fibra hay en el laminado (vs. resina). Para **laminado manual al vacío** el rango típico y alcanzable es:

| Proceso | Vf típico |
|---|---|
| Laminado manual sin vacío | 0.35 – 0.45 |
| **Laminado manual + bolsa de vacío** | **0.45 – 0.55** |
| Infusión de resina / autoclave | 0.55 – 0.65 |

**Usar Vf = 0.50 como valor de diseño.** Si el proceso entrega consistentemente un Vf más bajo (p. ej. 0.45 por exceso de resina), correr el script con `--vf 0.45` y verificar que el FSF siga siendo ≥ 1.50.

---

## 3. Cómo correr los scripts

### Requisitos

- Python 3.8 o superior
- Sin dependencias externas (solo librería estándar de Python)

### Paso 1 — Calcular propiedades del laminado

```bash
python3 inplaneG_v5.py --thickness 6.0 --beta 20 --vf 0.50 --json aleta_6mm.json
```

| Argumento | Valor | Significado |
|---|---|---|
| `--thickness 6.0` | 6.0 mm | Escala todas las capas al espesor objetivo |
| `--beta 20` | 20° | Rotación de orientación recomendada (ver [Glosario](#7-glosario)) |
| `--vf 0.50` | 0.50 | Fracción volumétrica de fibra |
| `--json aleta_6mm.json` | archivo | Guarda resultados para el siguiente script |

### Paso 2 — Verificar resistencia al flameo

```bash
python3 flutterEstimate_v3.py --json aleta_6mm.json
```

Eso es todo. El script lee el archivo generado en el paso 1 y produce el reporte de flameo.

### Opción: verificar con diferentes Vf

Si quieren evaluar el efecto de distintos Vf sin cambiar el apilado:

```bash
# Vf optimista (proceso bien controlado)
python3 inplaneG_v5.py --thickness 6.0 --beta 20 --vf 0.55 --json aleta_vf55.json
python3 flutterEstimate_v3.py --json aleta_vf55.json

# Vf conservador (proceso con algo de exceso de resina)
python3 inplaneG_v5.py --thickness 6.0 --beta 20 --vf 0.45 --json aleta_vf45.json
python3 flutterEstimate_v3.py --json aleta_vf45.json
```

---

## 4. Cómo leer los resultados

### Salida del script 1 (`inplaneG_v5.py`)

Al correrlo sin `--quiet`, imprime una tabla de barrido de ángulo beta. La línea más importante para ustedes es la de **beta = 20.0** (diseño recomendado):

```
  beta    D11       D22       D66       D16        G_xy     t      Nota
  [deg]  [N.m]     [N.m]     [N.m]     [N.m]     [GPa]   [mm]
  ──────────────────────────────────────────────────────────────────────
  20.0   ####.##   ####.##   ####.##   -##.####   #.####   6.000  RECOMMENDED — washout
```

**Columnas clave:**

| Columna | Qué mide | Para qué sirve |
|---|---|---|
| **D66 [N·m]** | Rigidez torsional | Entrada principal al cálculo de flameo — cuanto mayor, mejor |
| **D16 [N·m]** | Acoplamiento flexión-torsión | Debe ser **negativo** para beta=20° (efecto washout estabilizador) |
| **t [mm]** | Espesor total del laminado | Debe ser 6.000 mm al usar `--thickness 6.0` |
| **G_eff [GPa]** | Módulo torsional efectivo | Calculado como 12·D66/t³; entra directamente en la ecuación de flameo |

### Salida del script 2 (`flutterEstimate_v3.py`)

El reporte completo tiene esta estructura:

```
════════════════════════════════════════════════════════════════════════
  FIN FLUTTER ANALYSIS  [flutterEstimate_v3.py]
════════════════════════════════════════════════════════════════════════
  D66  = ####.#### N.m     G_eff = ##.#### GPa
  t    = 6.0000 mm

  FIN GEOMETRY
    cr = 300.0 mm   ct = 150.0 mm   b = 160.0 mm
    Λ  = 57.4 deg   AR = 0.7111   MAC = 233.33 mm
    t/c = 0.0257  (2.57%)

  ISA ATMOSPHERE  h = 1462.0 m
    ρ = 1.1049 kg/m³   T = 278.65 K   a = 335.94 m/s

  ROCKET  Vr = 652.7 m/s   M = 1.9420

  FLUTTER BOUNDARY CHAIN
    [1] NACA TN 4197 (subsónico)      : Vf = ####.# m/s   Mf = #.###
    [2] Corrección supersónica Ackeret : ×#.####  → Vf = ####.# m/s
    [3] Corrección por flecha (Λ=57.4°): ×#.####  → Vf = ####.# m/s

  FLUTTER SAFETY FACTOR
    V_flutter = ####.# m/s   V_rocket = 652.7 m/s
    FSF = #.####   [PASS ✓]   required ≥ 1.50   margin = +##.#%
```

### La línea más importante

```
    FSF = #.####   [PASS ✓]   required ≥ 1.50   margin = +##.#%
```

| Resultado | Interpretación | Acción |
|---|---|---|
| `[PASS ✓]` y FSF ≥ 1.50 | La aleta es segura aeroelásticamente | Proceder con la fabricación |
| `[FAIL ✗]` y FSF < 1.50 | La aleta puede flamear en vuelo | **No fabricar** — consultar al equipo de diseño |

**El "margin" indica cuánto margen adicional hay.** Un margin de +43% significa que la velocidad de flameo es 43% más alta que la velocidad del cohete en el punto de máxima presión dinámica.

### Cadena de correcciones de flameo

El script calcula la velocidad de flameo en tres pasos acumulativos:

| Paso | Método | Qué agrega |
|---|---|---|
| [1] NACA TN 4197 | Base subsónica | Estimación inicial conservadora |
| [2] Ackeret supersónico | Corrección M > 1 | Aumenta Vf porque en régimen supersónico la aerodinámica es menos acoplante |
| [3] Flecha (Λ = 57.4°) | Estabilización geométrica | La flecha trasera sube aún más Vf (factor ≈ 1.85) |

El valor final de **Vf después del paso [3]** es el que se compara con la velocidad del cohete para calcular el FSF.

---

## 5. Flujo completo de trabajo

```
┌─────────────────────────────────────────────────────────┐
│  DISEÑO Y VERIFICACIÓN — Aletas 6 mm Aconcagua          │
└─────────────────────────────────────────────────────────┘

  1. Definir parámetros de proceso
     └─ ¿Qué Vf espero lograr con mi bolsa de vacío? → 0.50

  2. Calcular propiedades del laminado
     └─ python3 inplaneG_v5.py --thickness 6.0 --beta 20 --vf 0.50 --json aleta_6mm.json
        └─ Verificar: t = 6.000 mm, D66 > 0, D16 < 0

  3. Verificar flameo
     └─ python3 flutterEstimate_v3.py --json aleta_6mm.json
        └─ Verificar: FSF ≥ 1.50  →  [PASS ✓]

  4. Si el proceso entrega Vf diferente al nominal
     └─ Repetir pasos 2 y 3 con el Vf medido (p.ej. --vf 0.47)
        └─ ¿Sigue siendo PASS? → OK para fabricar

  5. FABRICAR con el apilado de 22 capas definido en §2
     └─ Controlar espesor final con calibre → objetivo 6.0 ± 0.2 mm
```

---

## 6. Parámetros que pueden ajustar

Estos son los únicos parámetros que el equipo de manufactura necesita modificar según el proceso real:

| Parámetro | Comando | Rango útil | Cuándo cambiarlo |
|---|---|---|---|
| Fracción de fibra | `--vf 0.50` | 0.40 – 0.55 | Si el proceso entrega un Vf distinto al nominal |
| Espesor objetivo | `--thickness 6.0` | 5.5 – 6.5 | Si el herramental cambia el espesor final |

**No cambiar** `--beta`, `--layup`, `--cr`, `--ct`, `--span` ni `--sweep` sin consultar al equipo de diseño estructural.

---

## 7. Glosario

| Término | Definición |
|---|---|
| **Vf (fracción volumétrica de fibra)** | Proporción de fibra vs. resina en el laminado. Más fibra = más rígido. Para laminado manual al vacío: 0.45–0.55. |
| **D66 [N·m]** | Rigidez torsional del laminado. Resiste el giro de la aleta. Es el parámetro estructural que más influye en el flameo. |
| **D11, D22 [N·m]** | Rigideces a la flexión (en dirección del largo de la envergadura y de la cuerda). |
| **D16 [N·m]** | Acoplamiento entre flexión y torsión. Negativo en beta=20° → la aleta se tuerce en la dirección que descarga la carga aerodinámica (efecto washout, estabilizador). |
| **G_eff [GPa]** | Módulo torsional efectivo = 12·D66/t³. Entra directamente en la ecuación de velocidad de flameo. |
| **Beta (β)** | Rotación global de orientación de las fibras. Beta=0° es la orientación estándar. Beta=20° es el ángulo de adaptación aeroelástica recomendado para la geometría de Aconcagua. |
| **FSF (Factor de Seguridad al Flameo)** | Cociente entre la velocidad de flameo y la velocidad máxima del cohete. Debe ser ≥ 1.50. |
| **Flameo (flutter)** | Vibración aeroelástica autosustentada. Si ocurre, la aleta se destruye en milisegundos. |
| **Washout** | Efecto en el que la aleta, al doblarse hacia arriba, se tuerce reduciendo su ángulo de ataque. Estabilizador. Lo produce D16 < 0. |
| **Max-q** | Punto de máxima presión dinámica durante el vuelo (h ≈ 1462 m, M ≈ 1.94). Es el instante más crítico para el flameo. |
| **NCF biaxial** | Tela de fibra no-crimp (sin ondulación) en orientación ±45°. Máxima resistencia torsional. |
| **Tejido 0/90** | Tela tejida (woven) con fibras en 0° y 90°. Aporta rigidez a la flexión en ambas direcciones. |
| **ISA** | Atmósfera Estándar Internacional (ICAO). Modelo que calcula densidad, temperatura y velocidad del sonido según la altitud. |
| **Ackeret** | Corrección aerodinámica para régimen supersónico (M > 1). Eleva la velocidad de flameo estimada respecto al modelo subsónico. |
| **Flecha (Λ)** | Ángulo de barrido del borde de ataque de la aleta. En Aconcagua: Λ = 57.4°. La geometría de aleta con flecha trasera es inherentemente más estable al flameo. |

---

## Notas de fabricación — Laminado manual al vacío

- **Controlar el Vf durante el proceso:** el exceso de resina baja el Vf y reduce la rigidez. Usar rodillo desaireador firmemente. La bolsa de vacío debe mantenerse ≥ 0.85 bar durante todo el curado.
- **Respetar el orden de apilado:** el laminado es simétrico. Invertir el orden rompe la simetría y puede generar deformación por temperatura (warping) durante el curado.
- **Verificar espesor final:** medir con calibre en al menos 5 puntos. El objetivo es 6.0 ± 0.2 mm. Si el espesor promedio medido difiere, introducirlo en el script con `--thickness <valor_real>` y re-verificar el FSF.
- **La orientación ±45° de la tela DB300 está referenciada al eje longitudinal de la aleta** (dirección de la envergadura). Verificar marcas de referencia en el herramental antes de cortar.
