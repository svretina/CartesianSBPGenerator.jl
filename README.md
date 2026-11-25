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

`closure_depth` sets the number of boundary rows and `norm_bandwidth`
limits the largest allowed off-diagonal distance in the restricted block.
At the default depths of 7 and 9, respectively, the restricted-full
SBP¹(6,5) and SBP¹(8,7) trials give consistent systems with four and five free
parameters. The current norm selection heuristic leaves positivity for
these trials unresolved. Each design call tests a single closure configuration.

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
The current SBP generator still uses the direct closure equations. The
reduced norm-only construction is the next milestone.
