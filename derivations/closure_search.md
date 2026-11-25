# Closure and positivity search

For a fixed stencil radius, search increasing closure depths within each
allowed norm family. Explicit closure depths and sampled basis arrays fix
the closure size. Automatic norm selection tries the diagonal family first
when the theoretical order restriction permits it.

For a diagonal norm, first try the affine particular solution and a
candidate near the identity. In exact arithmetic, any remaining positivity
question reduces to strict linear inequalities in the free norm directions.
Fourier–Motzkin elimination combines lower and upper bounds to remove one
parameter at a time. Back-substitution picks a rational point inside each
open interval. A contradiction certifies failure for that configuration.
A work limit bounds elimination growth and leaves feasibility unresolved
if reached.

For a restricted full norm, check candidates with rational LDLᵀ in exact
arithmetic. Test a bounded set of displacements along independent norm
directions if the initial candidates fail. This finite search does not prove
infeasibility. Exact invariant nonpositive diagonal entries and a fixed
indefinite norm can certify infeasibility.

The diagonal regression targets are depths 11 and 14 for SBP¹(10,5) and
SBP¹(12,6), respectively. These match Table 1 of
[Albin and Klarmann](https://arxiv.org/abs/1403.5750v3).

The exact search reproduced both targets. It also reached a positive
restricted-full SBP¹(6,5) closure at depth 7. For restricted-full SBP¹(8,7),
depths 4–8 were inconsistent and depths 9–12 had no positive candidate
under the bounded norm selection procedure. The equal-order trials
SBP¹(2,2), SBP¹(4,4), and SBP¹(6,6) were inconsistent from their minimum
radii through depth 32 with minimum-radius stencils.
