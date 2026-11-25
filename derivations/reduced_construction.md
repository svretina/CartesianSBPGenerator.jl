# Reduced construction

The generator forms the local left coupling `C` from the centered stencil.
In dimensionless coordinates the closure satisfies `A*X = H*Y + K`, where
`A` is skew and `K = e₀e₀'X/2 - C*X_C`. The direct backend solves this
equation for the norm and skew entries together.

The reduced backend first solves for the norm entries using
`X'*(H*Y + K) + (H*Y + K)'*X = 0`. If the columns of `X` are dependent,
it also imposes `(H*Y + K)*Z = 0` for a basis `Z` of the kernel of `X`.
The norm is tested for positivity before the skew block is reconstructed.

The linear reconstruction expands the skew block in its independent entries.
The analytical reconstruction selects independent columns of `X`, evaluates
the pseudoinverse formula on those columns, and spans the remaining skew
freedom using antisymmetric outer products of vectors in the kernel of `X'`.
For rank `r` and depth `m`, that freedom has dimension `(m-r)(m-r-1)/2`.

Custom basis samples add constraints to the polynomial baseline. Their rows
cover the closure stencil, at the dimensionless nodes `0:closure_width-1`.
Reflection carries the same conditions to the right closure in its inward
coordinate. The interior stencil keeps its polynomial accuracy contract.

Regression tests compare the direct and reduced systems on the four
classical diagonal cases and the restricted-full SBP¹(4,3) case. They
check affine dimensions and exact membership in the other system, including
perturbations along the free skew directions. Another test gives a dependent
sample column for which symmetric compatibility alone is insufficient.
The numerical test reproduces the derivative of `exp(-ξ)` at each closure
in its inward coordinate and verifies the assembled SBP identity.

Automatic closure search, rank-adaptive enrichment, and operator
certification have their own milestones. In particular, the current
positivity candidate selection can leave a consistent system unresolved.
