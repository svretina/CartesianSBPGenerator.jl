# Certification of the period-two candidates

Run `julia --project=. experiments/certify_alternating.jl --pde` to reproduce
the campaign. Add `--dispersion-table` to emit both branches, phase and group
velocities, and smooth-mode overlap at every sampled wavenumber. Exact
[coefficients](alternating_coefficients.md) and the
[cut, parity, and spherical bulk derivations](alternating_compatibility.md)
are available separately.

## Locality and conditioning

All diagnostics below refer to the specific affine solutions selected by
the constructors. For conventional operators, D means interior order
`2p` with a diagonal norm and R means interior order `p+2` with a restricted-full norm.
All three have boundary order `p`.

| Boundary order | Family | Interior radius | Closure depth | Max boundary width | Nonzeros per row | Diagonal range of H/h | Cond(H) | Max abs(hD) |
| ---: | --- | ---: | ---: | ---: | --- | --- | ---: | ---: |
| 2 | Alternating | 2 | 5 | 7 | 4–5 | 0.389–1.096 | 2.82 | 5.27 |
| 2 | D | 2 | 4 | 6 | 2–4 | 0.354–1.229 | 3.47 | 1.74 |
| 2 | R | 2 | 3 | 5 | 3–5 | 0.333–1.375 | 4.25 | 2.00 |
| 4 | Alternating | 4 | 9 | 13 | 6–9 | 0.306–1.441 | 4.72 | 1.75e3 |
| 4 | D | 4 | 8 | 12 | 7–9 | 0.257–1.798 | 6.98 | 1.00e3 |
| 4 | R | 3 | 6 | 9 | 6–9 | 0.247–2.350 | 15.97 | 3.70 |
| 6 | Alternating | 6 | 15 | 21 | 7–15 | 0.099–1.881 | 19.01 | 9.76e7 |
| 6 | D | 6 | 14 | 20 | 7–14 | 0.00286–2.346 | 819.93 | 1.17e9 |
| 6 | R | 4 | 12 | 16 | 8–16 | 0.258–1.965 | 30.28 | 1.58e6 |

Both alternating phases have (2p) nonzero interior taps. The ratios
(w_B/w_A) are approximately 1.059, 1.155, and 1.288. The large
derivative coefficients are concentrated in the selected closures.
The norm condition number alone does not capture this problem.

## Bloch analysis

The 2,049-point scan covers (0 < θ ≤ π/2). The exact weighted
skew identity was additionally checked coefficient-by-coefficient.
At (θ≈0.0998) the physical phase errors were about
(-3.29e-6), (-1.16e-11), and (1.6e-15). The last is at
Float64 resolution. The respective physical group errors were
(-1.64e-5), (-1.05e-10), and (8.2e-15).
The folded phase velocities were about (-1.66), (-2.66), and (-3.48).
The sampled dimensionless bulk spectral radii were 1.373, 1.741, and 1.962.
The squared folded overlap of a smooth sampled mode at the same `θ`
was approximately `5.0e-9`, `1.1e-12`, and `1.6e-16`.
The same-order centered bulk radii were 1.000, 1.372, and 1.586,
and the centered `p+2` radii were 1.372, 1.586, and 1.731.
At the edge of the reduced zone, the respective alternating branch
frequencies were `±1.334`, `±1.528`, and `±1.573`. The sampled
folded smooth-energy fractions peaked there at roughly 0.020%,
0.130%, and 0.398%. The complete scan records both branch velocities.

At zero reduced wavenumber, both the constant and checkerboard vectors
are exact periodic-bulk zero modes. Every extracted stencil is odd in
its offsets, so both its even and odd tap sums vanish separately.
There are no additional zeros for the three certified bulk symbols.
For their antisymmetric stencils, factor `2 sin(θ)` from the real
matrix `-iM(θ)` and write each term using
`sin(l θ)/sin(θ)=U_(l-1)(cos θ)`. The resulting determinant is a
rational polynomial of degree 2, 6 or 10 in `cos θ`. Exact Sturm
counts found no roots on `[0,1]`, including the endpoints.
Reproduce this independent check with
`rg '^\| [246] \| [AB] \|' derivations/alternating_coefficients.md | python experiments/bloch_zero_certificate.py`.
The finite-domain exact nullspace tests below give nullity one.

## Grid and derivative refinement

The same coefficients were checked on all seven odd grids
`17:2:29`, `25:2:37`, and `37:2:49` for orders 2, 4, and 6, respectively.
Additional grids were used for the derivative and PDE refinements.
Checks include exact polynomial order at every row, positive diagonal
weights, SBP, constant-only nullspace, and reflection on the signed
grid. The phase-aligned transfer rejects an even grid size.

Five-grid tests used `exp(x)` and `sin(2πx)` on `[0,1]` in 192-bit
arithmetic. The grid interval counts were `(16,24,32,40,48)`,
`(24,36,48,60,72)`, and `(36,54,72,90,108)`. The five regions are
even and odd bulk phases, left and right closures, and the first two
bulk rows adjacent to each closure. The latter overlap the phase
groups. Final observed rates in that order were:

| Order | exp(x) | sin(2πx) |
| ---: | --- | --- |
| 2 | 1.86, 1.88, 2.05, 1.96, 1.86 | 1.92, 2.04, 1.91, 1.91, 1.43 |
| 4 | 3.85, 3.86, 4.06, 3.94, 3.85 | 3.99, 3.98, 3.85, 3.85, 2.71 |
| 6 | 5.84, 5.85, 6.06, 5.94, 5.84 | 5.98, 6.00, 5.85, 5.85, 4.30 |

The transition maxima for sine are sensitive to which coordinates are
sampled by those moving two-row windows. Their exact row-wise orders
still equal the requested order.

## SAT benchmarks and time-step limits

Advection uses exact sine inflow over eight crossing times. The
first-order wave system is evolved in characteristics with homogeneous
left inflow and a reflecting right wall. A compact smooth rightgoing
pulse reflects and is compared against the exact returning
characteristic at time one. The script checks the homogeneous
semidiscrete energy inequality and reports energy snapshots, phase
errors, and the checkerboard projection of the advection error.

Spatial errors use the matrix exponential of the linear generator.
RK4 uses the same step across families on each grid. A half-step run
measures its discrepancy from the semidiscrete result. These
discrepancies are affected by roundoff at the extremely small steps in
order 6. The RK4 CFL values below are spectral estimates for the two SAT
generators, not fully discrete bounds on transient amplification.

| Order | Finest common grid | Family | Spectral CFL | Advection max error | Reflected-wave max error |
| ---: | ---: | --- | ---: | ---: | ---: |
| 2 | 49 | Alternating | 0.364 | 2.63e-3 | 8.48e-2 |
| 2 | 49 | D | 2.07 | 7.10e-4 | 7.80e-2 |
| 2 | 49 | R | 2.07 | 5.49e-4 | 8.05e-2 |
| 4 | 73 | Alternating | 1.40e-3 | 3.38e-5 | 2.63e-2 |
| 4 | 73 | D | 2.89e-3 | 5.68e-6 | 2.23e-2 |
| 4 | 73 | R | 1.79 | 2.09e-6 | 2.69e-2 |
| 6 | 121 | Alternating | 5.52e-8 | 1.54e-6 | 1.43e-2 |
| 6 | 121 | D | 1.01e-8 | 1.45e-5 | 1.24e-2 |
| 6 | 121 | R | 3.14e-6 | 8.94e-8 | 9.75e-3 |

Alternating advection error decreased from `6.27e-2 → 8.61e-3 →
2.63e-3` at order 2 and `3.77e-3 → 1.84e-4 → 3.38e-5` at order 4.
At order 6 it decreased from `6.43e-5 → 1.61e-6 → 1.54e-6`.
That final change is too small relative to temporal and floating-point
effects to infer an advection rate. At the finest grids the alternating
checkerboard projections were about `2.45e-4`, `6.77e-6` and `1.06e-8`.
At order 6 the very large closures preclude a practical explicit step.

## Production recommendation

Keep the selected period-two families experimental. The exact closure and
parity results support further work, but the folded zero mode, large
closure coefficients and resulting SAT stability restrictions make
these particular choices unsuitable for default use. Preserve the
verified algebraic construction and investigate alternative admissible
closure choices before proposing promotion. No novelty claim is made.

The [external comparison report](external_validation.md) records the
independent cases and remaining reference work.
