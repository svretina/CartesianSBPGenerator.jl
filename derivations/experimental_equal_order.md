# Experimental equal-order SBP constructions

Run the campaign with `julia --project=. experiments/equal_order.jl`.
The experimental module is opt-in. The local centered-stencil generator
continues to enforce its existing moment obstruction.

## Construction and certification

The stencil experiments group upper-triangular entries of the
dimensionless norm and weighted derivative under reflection and, when an
entry touches the interior, by offset and periodic phase. Boundary
weighted-derivative blocks can be dense within the selected closure.
The norm bandwidth can be set separately inside those blocks. Antisymmetry fixes the
lower part of the weighted derivative, with the two endpoint diagonal
entries coming from half of the boundary matrix. The exact linear system
imposes the polynomial equations at all nodes. For a repeating compact
norm, it also normalizes the interior norm row sum to one.

The search tests a bounded set of affine candidates for positivity.
For each accepted norm the code constructs the derivative by a solve and
checks the resulting row-wise exactness, SBP identity and nullity. If the
linear system is consistent and the candidate search fails, the result is
unresolved. The affine solution stays available for further investigation.

The wider centered systems leave the interior coefficients unknown and
fix the interior norm to identity. For orders 2, 4 and 6 and each tested
extra radius, the first omitted moment vanished on both the particular
solution and every affine direction. The exact Lagrange element cases
provide positive reference operators of those same orders at all nodes.

## Results

These fixed-grid candidates achieved the requested order at every row
and had exact nullity one. Each was also checked after transferring the
same coefficients to two finer grids.

| Order | Interior structure | Closure depth | Grid sizes | Norm condition number |
| --- | --- | ---: | --- | ---: |
| 2 | Explicit, period 2, diagonal norm | 5 | 17, 33, 65 | 2.82 |
| 4 | Explicit, period 2, diagonal norm | 9 | 25, 49, 97 | 4.72 |
| 6 | Explicit, period 2, diagonal norm | 15 | 37, 73, 145 | 19.01 |
| 4 | Compact, pentadiagonal norm | 6 | 17, 33, 65 | 8.32 |
| 6 | Compact, norm half-bandwidth 3 | 10 | 25, 49, 97 | 28.34 |

On the finest grid the nonzero counts of `H`, `Q` and `D` were
`(65, 264, 264)`, `(97, 774, 774)`, `(145, 1732, 1732)`,
`(319, 392, 4224)` and `(667, 806, 9408)`,
respectively. A direct application of the explicit periodic operator can
use linear work with sparse storage. Direct multiplication by the stored
compact derivative takes quadratic work. Applying the compact relation
with a reusable banded factorization would take linear work for fixed
bandwidth. The experiment stores dense matrices for exact certification,
and the campaign prints measured construction times.

For the manufactured function `exp(x)` on `[0,1]`, the final refinement
gave boundary rates about 1.96, 3.94 and 6.61 for the explicit cases and
3.96 and 6.07 for the compact cases. The sixth-order explicit interior
error fell to about `8.3e-14`. The corresponding floating-point rate is
sensitive to roundoff. These are observed rates over the tested grids.

The search also found diagonal-norm cases with period 4 at all three
orders, and with period 3 at order 2. Several smaller higher-order
explicit configurations had inconsistent equations, while others were
consistent with positivity unresolved.
The tested tridiagonal order-4 cases were inconsistent. For order 6,
some pentadiagonal systems were consistent with positivity unresolved.
The successful order-6 example uses a wider norm.

## Literature and software

- [Mattsson and Ljungberg Rydin (2023)](https://doi.org/10.1016/j.jcp.2022.111743)
  supplies the banded-norm compact framework. Its published implicit operators
  use their own boundary orders. This experiment solves different equal-order
  constraints.
- [Del Rey Fernández, Boom and Zingg (2014)](https://www.sciencedirect.com/science/article/pii/S002199911400076X)
  discusses generalized nodal SBP with nonrepeating interior stencils.
- The [authors' MATLAB repository](https://github.com/ylvarydin/imlicit_sbp)
  was inspected at commit `4219a72406f94356597c98490b9d4334644bf5e2`.
  Its first-derivative source `SBP_I6.m` uses banded interior matrices and
  separate closure blocks. The checkout has no license file. No code or
  coefficient tables were copied.
