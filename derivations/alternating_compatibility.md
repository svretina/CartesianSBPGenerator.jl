# Period-two compatibility and signed-radius parity

Write the dimensionless derivative coefficients as `d[r,l]` for
phase `r = A,B` and offset `l`. Phase A is the even zero-based index.
With `q[r,l] = w[r] d[r,l]`, the interior SBP condition is
`q[r,l] = -q[r+l,-l]`, where the phase index is taken modulo two.
The Taylor equations on each row are
`Σ_l d[r,l] l^k = δ_{k1}` for `0 ≤ k ≤ p`.

## Cut identity and leading error

For any cut between nodes `k` and `k+1`, apply polynomial exactness
to the rows to the left of the cut, paired with the coefficient matrix
of `(y-x)^p`. The diagonal norm terms vanish on the diagonal and the
interior skew terms cancel. The endpoint term vanishes because
`(x-x)^p=0`. The remaining sum is
`F_p(k) = Σ_{i≤k<j} Q[i,j] (j-i)^p = 0`.
For a cut in the repeating region, averaging `F_p` over both phases
counts each positive offset `l` exactly `l` times over two cuts,
with each phase represented equally. Antisymmetry of `Q`
and odd parity of the power `p+1` then give
`w_A μ_A + w_B μ_B = 0`. This is an exact necessary identity.

The selected candidates have both phase stencils antisymmetric.
Their positive odd-offset weighted coefficients coincide.
For even offsets `l=2,4,…,p`, set
`v_l=2l(q[A,l]-q[B,l])` and `z_l=l²`. The difference of the
Taylor equations gives `Σ_l v_l z_l^k=(w_A-w_B)δ_{k0}` for
`0 ≤ k < p/2`. Evaluating at zero the interpolant of
`z^(p/2)` on these nodes gives
`w_A μ_A - w_B μ_B = (w_A-w_B)(-1)^(p/2+1) Π_l l²`.
The weights differ, so both weighted moments are nonzero and opposite.
This argument explains the leading error in the constructed operators.
It does not establish positive-norm existence at arbitrary order.

## Fourier symbol and phase alignment

In the node gauge, the dimensionless symbol entry `M[r,t]` sums
`d[r,l] exp(i l θ)` over offsets taking phase `r` to phase `t`.
The cell-gauge symbol in the request is unitarily equivalent to this
one by the phase factor on odd nodes. Weighted skew symmetry gives
`W M + M* W = 0` for every real `θ`. The Hermitian matrix
`-i √W M (√W)^(-1)` yields real frequencies. Its differentiated
Rayleigh quotient gives group velocity on a simple branch.

Applied to the sampled smooth mode `(1,1)`, the first signed bulk
Taylor defect is proportional to `(C_A,C_B)`. Its weighted projection
onto `(1,1)` vanishes by the cut identity. Thus that first defect
couples to the folded direction instead of shifting the smooth
Rayleigh quotient directly. This is consistent with the small
low-wavenumber physical-branch errors in the numerical scan.

Each extracted stencil is odd in its offsets. For an origin on a grid
node, reflection preserves the phase while changing the sign of the
offset. Both bulk phases therefore satisfy origin parity. The exact
finite-domain constructions reflect about their central node on odd
grids. Seed coefficients are transferred only to odd grids with the
same reflection phases. Even grids reflect about a face and swap the
phases. They need a separately assembled closure or another origin
layout.

Adding node pairs only translates the row templates in dimensionless
coordinates. Once all phases and the two closures are represented
without overlap, translation preserves polynomial exactness, diagonal
weights remain positive, and the grouped entries preserve SBP.
These arguments transfer the local checks to larger phase-aligned grids.
The constant-only finite-grid nullspace was checked on the stated grids.
The local argument alone does not prove the rank for all grid sizes.

## Spherical bulk implication

Away from the origin and boundary, choosing both mass matrices as
`diag(h w_i r_i²)` gives the candidate divergence
`R⁻² G R²` by the Cartesian adjoint identity. For interior radius
`r_i > 0`, test row phase `a` on local monomials by expanding
`(r_i + j h)²` and using the moments of that phase. In particular,
monomial `r^k` is exact for `k+2 ≤ p`. For smooth data at a fixed
nonzero radius the leading error is
`h^p C_a (r² f)^(p+1)(r_i)/r_i²`. The factor `1/r_i²`
precludes carrying that fixed-radius estimate uniformly through the
origin. For smooth odd radial data the first potentially missed odd
power is `r^(p-1)`. At the first few points `r_i=O(h)` its error
can be `O(h^(p-2))`.

Preserve both phase weights in any bulk quadrature reuse and match the
phase at each closure. Origin regularization and its additional
polynomial conditions remain necessary. The reflected bulk stencil
itself needs no parity transition at a node-centered origin. This
derivation does not validate a spherical origin or scri closure.

I inspected the sibling implementations to check what can be reused.
In `SphericalSBPOperators.jl`, `src/core/fullgrid.jl` currently obtains
the signed Cartesian grid from an upstream derivative source.
The split-mass routines in `src/non_diagonal_mass/sbp4.jl` already
derive their fixed bulk entries from the projected Cartesian mass.
That algebra can retain its diagonal bulk structure and pick up both
phase weights from an alternating Cartesian mass. The constructor
needs an input path for the experimental gradient and mass, and the
closure row ranges need to be checked against the longer stencil.
In `SSBP.jl`, the mass-first system fixes a bulk scale of one relative
to the supplied volume weights. Reuse would require supplying the
phase-dependent baseline before forming the fixed moment contribution.

This does not justify reusing the closure equations with their old
row counts. For example, the split fourth-order solver imposes
`D*r = 3` away from the outer closure and `D*r^3 = 5r^2`
near the origin. Exact Cartesian order four suffices for the first
identity in the bulk, while the second requires exact differentiation
of `r^5` and is not supplied by the unchanged bulk formula.
The analogous split sixth-order test of `D*r^5` needs differentiation
of `r^7` near the origin. The second-order Cartesian case already
fails bulk exactness on `D*r = 3` for this uncorrected mass ansatz.
For the first missed odd monomial `f=r^(p-1)`, the exact bulk defect
is `μ_a h^p / r_i²`. No origin or outer mass solve was changed here.

The diagonal Cartesian weights integrate polynomial products through
degree `2p-1`, from the same discrete integration-by-parts identity.
On a symmetric signed grid, folding the weights gives the corresponding
half-interval quadrature on even data. After multiplication by `r²`,
the even radial monomials through degree `2p-4` are still integrated
exactly by that baseline. The entry at the origin is zero and must
be replaced as part of the regularized mass construction.
