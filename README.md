# CartesianSBPGenerator.jl

`CartesianSBPGenerator.jl` constructs finite-difference and, progressively,
summation-by-parts (SBP) operators on Cartesian grids. Derivative order,
interior accuracy order, and boundary accuracy order are independent inputs.

The current implementation completes the generic centered-stencil milestone,
provides the exact affine linear-algebra core, and constructs the established
diagonal-norm first-derivative families SBP¹(2,1), SBP¹(4,2), SBP¹(6,3), and
SBP¹(8,4) from their defining equations. It also constructs the restricted-full
SBP¹(4,3) family. Independent second derivatives remain under development.

```julia
using CartesianSBPGenerator

stencil = centered_stencil(
    derivative_order = 2,
    accuracy_order = 6,
    radius = 3,
    arithmetic = Rational{BigInt},
)

coefficients(stencil)
# 1//180 .* [2, -27, 270, -490, 270, -27, 2]
```

Wider stencils retain their free parameters instead of choosing hidden
coefficients:

```julia
family = centered_stencil(
    derivative_order = 1,
    accuracy_order = 4,
    radius = 3,
)

nullity(family)                     # 1
coefficients(family; parameters=[0])
```

Construct and instantiate a diagonal-norm SBP operator with exact rational
coefficients:

```julia
family = design_sbp(4, 2; derivative_order=1, norm=:diagonal)
op = instantiate(family; N=101, domain=(0, 1))

D = derivative_matrix(op)
H = mass_matrix(op)
B = boundary_matrix(op)

H * D + transpose(D) * H == B
```

Use `design_sbp(4, 3; norm=:restricted_full)` for a fourth-order interior
stencil with a third-order boundary closure. The solver recovers a
three-parameter affine family, matching the parameter count reported by
[Strand (1994)](https://doi.org/10.1006/jcph.1994.1005).
The selected member has an exactly certified positive-definite norm. Inspect
its dimensionless left block with `family.left_norm`.

`closure_depth=:auto` searches increasing closure depths up to
`max_closure_depth` (32 by default). An integer closure depth tests that
one size. `norm_bandwidth` limits the largest off-diagonal distance in a
restricted block. `positivity_margin` requires the norm minus that
margin times the identity to be positive definite.

The exact search reproduces the diagonal SBP¹(10,5) and SBP¹(12,6) closure
depths of 11 and 14 and finds restricted-full SBP¹(6,5) at depth 7.
`family.diagnostics.attempts` lists the tried configurations.
`SBPConstructionError.attempts` distinguishes inconsistent equations,
certified positivity infeasibility, and unresolved norm selection.
See the [closure search notes](derivations/closure_search.md).

For this local centered-stencil construction, a boundary order at least
as high as the exact interior order forces the interior stencil's first
omitted moment to vanish. The default reports this obstruction before
searching. To explicitly permit a higher interior order while retaining
the requested boundary order, use
`design_sbp(4, 4; norm=:restricted_full, order_policy=:allow_upgrade)`.
The resulting family records the original orders separately from the
constructed orders. See the [moment argument](derivations/equal_order_obstruction.md).

Run the tests with:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

See [the current specification](specification.md) for the revised
first-derivative milestone sequence and [the literature audit](derivations/literature.md)
for the construction and validation references.

The revised algebra foundation provides `PolynomialBasis`, `MonomialBasis`,
and `BasisSamples` for value and derivative samples. `solve_constraints`
uses exact rational elimination by default and accepts `backend=:numerical`
for a LinearSolve QR solve with independent SVD rank and residual checks.
The default `construction=:reduced` solves for the boundary norm before
reconstructing the skew part of the weighted derivative. Use
`construction=:direct` to solve for both together or
`skew_reconstruction=:analytical` to use the pseudoinverse formula.
See the [reduced construction notes](derivations/reduced_construction.md).

Additional accuracy functions can be supplied as value and derivative
samples at the dimensionless closure nodes. For example, the following
adds `exp(-ξ)` to the left closure of a second-order interior operator
with six boundary rows and a radius-one interior stencil:

```julia
ξ = 0:6
values = reshape(exp.(-ξ), :, 1)
basis = BasisSamples(values, -values, ["exp(-ξ)"])
family = design_sbp(2, 1; closure_depth=6, basis, arithmetic=:numerical)
```

The right closure uses the reflected inward coordinate. The additional
conditions apply to the closures, while the interior stencil retains its
prescribed polynomial accuracy.

## Construction methods

- `design_sbp(4, 2; method=:classical)` solves the norm and skew equations
  together. `construction=:direct` is the equivalent selector.
  Add `arithmetic=:numerical` to use the LinearSolve backend.
- `design_sbp(4, 2; method=:reduced)` solves the norm-only system first.
  Both stencil methods return resolution-independent families.
- `generate_sbp(basis; method=:function_space)` builds an operator on one
  grid and enforces exactness on the supplied basis at every node. Include
  physical derivatives in `BasisSamples`; the constructor adds constants.

For example, to differentiate a trigonometric space on `[0, 1]`:

```julia
x = range(0.0, 1.0; length=7)
s, c = sinpi.(2 .* x), cospi.(2 .* x)
basis = BasisSamples(hcat(s, c), 2π .* hcat(c, -s), ["sin", "cos"])
op = generate_sbp(basis; method=:function_space, arithmetic=:numerical, bandwidth=2)
D, H, Q = derivative_matrix(op), mass_matrix(op), weighted_derivative_matrix(op)
D * basis.X ≈ basis.Y
```

The fixed-grid method supports diagonal and restricted-full norms and
reflects both ends of the grid. `bandwidth` constrains the weighted
derivative. The default permits a dense weighted derivative.
See the [function-space derivation](derivations/function_space.md).

## Experimental equal-order construction

The opt-in `CartesianSBPGenerator.Experimental` module can construct
alternating diagonal-norm operators analytically from conventional SBP
families. The resulting derivative has exact order `p` at every row.
Use `p = 2, 4, 6` and an odd grid size with room for both closures and
a full interior period.

```julia
using CartesianSBPGenerator: Experimental
op = Experimental.analytical_perturbation(4; N=41)
op.diagnostics.bulk_taps
```

See the [analytical derivation and matched-grid comparison](derivations/analytical_perturbation.md)
for the construction, reproducible benchmarks, and stencil plots.

The experimental module also includes the original direct periodic
construction, compact banded norms, wider centered stencils with unknown
coefficients, and exact polynomial elements.

```julia
using CartesianSBPGenerator: Experimental
system = Experimental.assemble_stencil_experiment(
    4; N=25, closure_depth=9, stencil_radius=4, period=2,
)
result = Experimental.search_stencil(system)
refined = Experimental.refine_stencil(system, result.operator; N=49)
```

This period-two example has a diagonal norm and a sparse explicit
derivative of exact fourth order at every row, including the boundary.
See the [experiment report](derivations/experimental_equal_order.md) for
the sixth-order result, compact alternatives, and refinement checks.

The [alternating certification report](derivations/alternating_certification.md)
includes exact interior coefficients, Bloch dispersion and SAT benchmarks
for the original direct construction. Its high-order boundary closures
impose severe explicit time-step restrictions.
