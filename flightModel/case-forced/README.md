# case-forced/ — case template, forced motion of a fin

`rhoPimpleFoam` in physical time on a moving mesh: fin 0 moves as a sum of the
structural modes in `constant/flutterModes` times an imposed signal, and
`genForces` writes the generalised aerodynamic forces at each time step. Parent
`case-pimple/` (file `parent`); `pimple.exclude` lists the files of the parent
it drops.

> **CAUTION** Not usable yet: on h0 the time-accurate run diverges after 6 to 8
> time steps, also with a static mesh. `fin_flutter/docs/PROCEDURE.md`, step 7.

```bash
../newCase.sh ~/runs/f12-h0-p1 --mesh ~/meshes/h0 --solver forced \
    --Minf 1.2 --pInf 60547 --Tinf 280.15 \
    --init ~/runs/m12-a00-h0-pimple --flutterModes ~/flutter/aconcagua/pulse-m1-A0.5mm
```

The procedure, the files and the post-processing: `../../fin_flutter/docs/PROCEDURE.md`.
