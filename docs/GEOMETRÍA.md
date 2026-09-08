`mesh/aconcaguaGeom.py` es la **geometría de registro** del Aconcagua. 
Prescindimos del .stl para lograr mejor calidad de elementos y tener una malla estructurada.
## Parámetros base

```python
L_NOSE   = 0.80000   # m  longitud del cono (Von Kármán / LV-Haack, C = 0)
L_CYL    = 2.03000   # m  cilindro
L_TAIL   = 0.12500   # m  boattail (tronco de cono recto)
R_BODY   = 0.07550   # m  radio del cilindro  (D = 0.151)
R_BASE   = 0.05500   # m  radio de la base
```

## Geometría de aletas

```python
N_FINS       = 4          # fijado por la topología de la malla, no se cambia
FIN_T        = 0.012      # m  espesor total
FIN_ROOT_R   = 0.07500    # m  radio de la cuerda raíz (enterrada 0.5 mm en el cuerpo)
FIN_TIP_R    = 0.23550    # m  semi-envergadura
FIN_ROOT_LE  = X_BASE - 0.42578
FIN_ROOT_TE  = X_BASE - 0.12531    # = X_BODY_2: la raíz termina donde empieza el boattail
FIN_TIP_LE   = X_BASE - 0.17500
FIN_TIP_TE   = X_BASE - 0.02500    # la punta termina sobre el boattail
FIN_LE_BEVEL = 0.03320    # m  bisel del borde de ataque (ABSOLUTO, no % de cuerda)
FIN_TE_BEVEL = 0.02116    # m  bisel del borde de fuga
```
## Cambiar la geometría

1. Editar las constantes en `aconcaguaGeom.py`.
2. Correr archivo como `python aconcaguaGeom.py`. La función `validate()` chequea el cierre del cono sobre el cilindro, monotonía de la ogiva (el mapa del casquete debe ser biyectivo),
   raíz de la aleta dentro del cuerpo, biseles que no se solapan en la punta, `N_FINS == 4`. Estos chequeos son **editables**.
	1. Si hay un STL, `python aconcaguaGeom.py ruta.stl` compara `r_body` contra sus vértices.
3. `python build.py --plan`: `validate_params()` de la malla vuelve a chequear
   los acoples con la geometría (`ZONE_R[0] > R_BODY`, `FIN_X_LEAD` sobre el cilindro, etc.) y nombra lo que hay que ajustar.

### Qué parámetros de malla suelen moverse con la geometría?

| Cambio                      | Revisar                                                                                                |
| --------------------------- | ------------------------------------------------------------------------------------------------------ |
| `R_BODY` mayor              | `ZONE_R[0]` (tiene que quedar fuera del cuerpo), `ZONE0_R_WAKE`, `FIN_TIP_R`                           |
| cuerpo más largo            | `UPSTREAM_L`/`DOWNSTREAM_L` son en longitudes de cuerpo: el dominio crece solo. `X_WAKE_2` es absoluto |
| aletas más largas en cuerda | `FIN_X_LEAD` (el bloque refinado tiene que empezar sobre el cilindro)                                  |
| aletas más altas            | una zona radial fuera de `FIN_TIP_R` (`--preset fintip` como modelo)                                   |
| otro perfil de cono         | reemplazar `r_nose`; mantener monotonía y `r_nose(L_NOSE) == R_BODY`                                   |

## Valores de referencia

`reference_values(fraction)`: `Aref = π D²/4 × fraction` (sección del cuerpo
por la fracción de 360° que hay en la malla), `lRef = D`. Convención de
cohetería (Barrowman, OpenRocket).
