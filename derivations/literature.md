# Mathematical references

## Local compatibility and existence

Albin and Klarmann §2, Theorems 1 and 2, eliminate the skew boundary unknowns
before solving for the norm. Their `(s,t,r)` maps to interior order `2s`,
boundary order `t`, and closure depth `r`. The interior stencil determines the
boundary-to-interior coupling. Equation (9) is the block-norm compatibility
system, and equation (10) is its diagonal moment specialization.
[Albin and Klarmann, arXiv v3](https://arxiv.org/abs/1403.5750v3)

For a closure basis matrix `X` and target `F`, the equation `A*X = F` with
`A' = -A` requires both `X'*F + F'*X = 0` and `F*ker(X) = 0`.
For sufficiency, restrict to a basis of the column space of `X`. The kernel
condition makes the map from the columns of `X` to those of `F` well-defined.
Skew compatibility makes its restriction to that column space skew, and the
transpose block fixes its coupling to the orthogonal complement. The skew
block on the complement remains free. This covers the one-row SBP¹(2,1) case.

Section 3 separates exact equality feasibility from positivity. A unique
nonpositive norm rules out that configuration. An indefinite particular
solution in an affine family calls for a feasibility search. Table 1 gives
diagonal reference closure depths 1, 4, 6, 8, 11, and 14 for interior orders
2, 4, 6, 8, 10, and 12 under the paper's matching boundary orders and ansatz.

## Sampled bases

Glaubitz, Nordström and Öffner §4 use sampled values and derivatives in
`Q*X = H*Y`. Their equation (4.4) gives the symmetric compatibility identity.
Theorems 4.4 and 4.5 connect diagonal norms to quadrature exact on derivatives
of basis products, with a full-column-rank hypothesis for sufficiency.
The local Cartesian system includes the known coupling across the closure cut.
[Glaubitz, Nordström and Öffner, arXiv v2](https://arxiv.org/abs/2203.05479v2)

## Operator quality

The 2026 preprint's v2 title is *Beyond the summation-by-parts property:
nullspace consistency, sparsity, and regularization for FSBP operators*.
Proposition 4.4 explains how an unrestricted dense minimum-norm selection can
have rank at most twice the basis dimension when constants are included.
Section 4 links additional null modes to stationary errors and identifies
poor resolution of physical modes as a separate mechanism.
[Glaubitz et al., arXiv v2](https://arxiv.org/abs/2602.10786v2)

For this package, check exact rank of the weighted derivative on odd and even
grids, measure small singular values, and report the eigenvalues of the
inflow-SAT derivative. The endpoint-decoupled mass block makes the inflow
penalty a single diagonal correction. Compare observed convergence by region
and validate the energy identity independently.
