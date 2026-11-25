# Independent reference comparisons

The normal test suite uses frozen fixtures and has no dependency on the
reference repositories or external construction tools.

- The pinned upstream `SummationByPartsOperators.jl` v0.5.106-DEV
  generates exact Mattsson–Nordström diagonal closures at interior orders
  2, 4, 6, and 8. The test checks norm structure, polynomial exactness,
  the boundary identity, and interior taps independently. The order-2
  and order-4 closures match both matrices exactly. The higher-order
  closures satisfy our exact direct and reduced systems. Solving for
  their free skew coordinates reconstructs the upstream coefficients
  exactly with its norm fixed.
- Two upstream dense function-space fixtures use a fixed seven-node
  uniform grid with linear or trigonometric bases. The source solver
  uses Optim and ForwardDiff only during fixture generation.
  This package verifies the external norm and basis conditions, then
  compares its algebraic reconstruction at the same fixed norm via the
  skew null directions. The tests reproduce the singular values and SAT
  spectra and distinguish the linear fixture with four null directions
  from the trigonometric fixture with one at a relative rank tolerance
  of `1e-10`. The original diagnostic function from the 2026 checkout
  independently reproduced these ranks and their tolerance sensitivity.
  Re-run it with
  `julia experiments/external_diagnostics.jl <2026_SBP_not_enough checkout>`.
- Exact closure-search comparisons against Albin–Klarmann's published
  diagonal minimum depths are documented in [the search report](closure_search.md).
- The original MATLAB `compute_QF`, `generate_span`, and `compute_FSBP`
  also ran in GNU Octave 11.3.0. The polynomial case uses three uniform
  Lobatto nodes and the quadratic basis. It compares the derivative and
  mass against the unique exact solution from our function-space
  constructor. The trigonometric case uses four uniform nodes and the
  original trapezoidal mass. Both external skew reconstructions match
  ours at the same fixed mass within numerical tolerance.

Regenerate the fixtures by activating a temporary Julia environment and
developing the checked-out upstream package at
`20f6b1293622fd6afef062b5fc8780b5b617af47`. Run
`experiments/external_classical.jl` and
`experiments/external_function_space.jl` and capture their standard output.
For the second script install Optim 2.3.2 and ForwardDiff 1.5.0 in that
environment. The frozen data include both generated derivative matrices
and their norm information, in compressed closure format for the
repeating classical stencils.

For the MATLAB fixture, set `FSBP_REFERENCE` to the authors' checkout
at `e82d787d7c61b89d7aa98715cbd294d726ec28b8`. Put `experiments`
on the MATLAB or Octave path, then run `external_matlab_fsbp` and
capture standard output. Octave used the local pseudoinverse shim
for the missing MATLAB `lsqminnorm` function. The Octave container
manifest was
`sha256:185db7993e000d4f3f6e7bbbf7fb3f999f52e799ea52231ad8a15353381e0dcb`.

The broader external campaign remains open. This pass did not run
the full 2026 reproduction campaign. Only the 2026 diagnostic routine
was executed on the frozen fixtures.
It did not find an independently matched restricted-full coefficient
fixture, and it did not validate the future automatic additional-accuracy
policy against external data. The independent classical and dense FSBP
fixtures do not directly certify the new alternating families.
