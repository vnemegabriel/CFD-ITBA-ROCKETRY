# Supersonic Fin Flutter Solver

Aeroelastic stability of a composite CFRP rocket fin (IREC 2026, Team 207, ITBA): CLT laminate → MITC4 Mindlin plate FE → modal basis → supersonic quasi-steady aerodynamics → exact state-space stability over the flight envelope.

## Usage

```bash
pip install -e .[test]
python -m finflutter cases/ar1_8mm.yaml   # prints summary, writes results/ar1_8mm.json (tagged with git commit)
pytest tests
```

The run is configured entirely by the case YAML: geometry, mesh, laminate (`layup`, `thickness_mm`, `vf`, `beta_deg`), clamped root intervals `clamp_x`, mode count, loss factor `g`, aero theory (`ackeret` | `lighthill`, `sides`), and the flight CSV.

## Layout

| Module | Content |
|---|---|
| `finflutter/laminate.py` | Halpin-Tsai + CLT for DB300/GA90R layups, closed-mould Vf, β tailoring |
| `finflutter/structure.py` | Fin mesh, MITC4 plate K/M, aero operators ∫NᵢNⱼ and ∫Nᵢ∂Nⱼ/∂x, BCs, modes |
| `finflutter/aero.py` | GAF providers `(mach) → (A0, B1)`: `PistonTheory`; `TabulatedGAF` for CFD/panel-code Q(M, k) |
| `finflutter/aeroelastic.py` | `η̈ + (gΩ + qB1/U)η̇ + (Ω² + qA0)η = 0`; critical q by sweep + bisection; flutter/divergence classification |
| `finflutter/flight.py` | ISA, flight CSV (time, altitude ft, Vz m/s) |
| `tests/test_physics.py` | Leissa CFFF frequencies, Dowell SSSS panel flutter (λ ≈ 512), aero-damping sign, GAF fit, CLT regression |

## Model

Pressure jump across the fin (flow +x, both faces): Δp = sides · c₀ q (∂w/∂x + (c₁/c₀) ẇ/U), with
`ackeret`: c₀ = 2/β, c₁ = c₀(M²−2)/(M²−1); `lighthill`: c₀ = c₁ = 2/M. Generalised force −∫φᵢΔp dA.

Margins are reported at constant Mach as q_crit/q_flight (velocity margin = √ of it).

## Limitations

- Below M = √2 the Ackeret damping coefficient is negative, and it predicts single-mode flutter there; Lighthill piston theory needs M² ≫ 1. Neither is valid in the M 1.05–1.4 band, which needs unsteady linear supersonic theory or CFD-derived GAFs (`TabulatedGAF`).
- No thickness (2nd-order piston) term, no fin–body interference, structural damping g assumed.

## References

| Source | Used for |
|---|---|
| Lighthill (1953) *J. Aero. Sci.* 20(6) | Piston theory |
| Dowell (1970) *AIAA J.* 8(3); Dowell, *Aeroelasticity of Plates and Shells* | Supersonic quasi-steady aero, panel flutter benchmark |
| Leissa (1969) NASA SP-160 | Cantilever plate frequencies |
| Bathe & Dvorkin (1985) *IJNME* 21 | MITC4 element |
| Weisshaar (1981) *J. Aircraft* 18(8) | Aeroelastic tailoring |
| Jones (1999) *Mechanics of Composite Materials*; Halpin & Tsai (1969) | CLT, micromechanics |

## License

MIT. See [LICENSE](LICENSE).
