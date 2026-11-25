# External implementation audit

| Repository | Inspected commit | License |
| --- | --- | --- |
| [FSBP](https://github.com/jglaubitz/FSBP) | `e82d787d7c61b89d7aa98715cbd294d726ec28b8` | BSD-2-Clause |
| [2026_SBP_not_enough](https://github.com/JoshuaLampert/2026_SBP_not_enough) | `1cd3b07d6012e9bf88bb54e02289af1986c0f35a` | MIT |
| [SummationByPartsOperatorsExtra.jl](https://github.com/JoshuaLampert/SummationByPartsOperatorsExtra.jl) | `3aa493678435cf76aa67ca7cf7b28e0e98bfbcea` | MIT |
| [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl) | `20f6b1293622fd6afef062b5fc8780b5b617af47` | MIT |
| [LinearSolve.jl](https://github.com/SciML/LinearSolve.jl) | `1c05d94cf55ded0e5bdccdea93594675646877b8` | MIT |

- FSBP's `generate_span.m` dispatches to basis-specific evaluators.
  `generate_span_poly.m` samples monomials and derivatives and supplies
  product-derivative moments for equation (4.6). `compute_QF.m` selects
  quadrature and `compute_LSQF.m` grows the node set until the least-squares
  weights are positive within tolerance. `compute_FSBP.m` vectorizes
  `S*X = H*Y - B*X/2` and appends skew constraints. Use the sampled-basis
  interface and the skew equation as references. `solve_linear_SAT.m`
  provides a periodic advection and energy-monitoring example.
- The 2026 experiments' `code/operator_diagnostics.jl` independently computes
  weighted and unweighted residuals, compares SVD rank thresholds, and applies
  a unit inflow penalty before computing SAT eigenvalues.
  `code/stationary_modes.jl` projects advection errors onto stationary modes.
  `code/dispersion_analysis.jl` uses a Bloch matrix to separate physical
  dispersion and stationary modes. These are validation references.
- Extra.jl's `src/function_space_operators.jl` accepts functions and sampled
  nodes and dispatches to backend-specific construction.
  `src/utils/moments.jl` builds boundary product moments.
  `test/test_FSBP_reproducing.jl` compares both derivative and mass matrices
  against independently generated operators. Use these interface and
  regression-test patterns. Its construction backends depend on optimization.
- SummationByPartsOperators.jl's `src/SBP_operators.jl` stores local
  coefficients with a grid and spacing factor and applies them via `mul!`.
  Use its operator interface and established coefficients as independent
  regression references. [Issue #130](https://github.com/ranocha/SummationByPartsOperators.jl/issues/130)
  links to the Albin–Klarmann paper and includes no standalone implementation.
- LinearSolve.jl's `src/factorization.jl` provides pivoted QR and SVD
  algorithms. The numerical backend can use a rank-revealing factorization
  for the affine basis, a LinearSolve particular solve, and an independent
  residual check. The rational backend remains independent.

The repositories were cloned outside the package tree. The implementation
uses independently written algorithms.
