# Fixed-grid function-space method

The classical and reduced methods construct reusable local stencils for
prescribed polynomial orders. The function-space method uses fixed samples
of values and physical derivatives on a uniform grid. It appends the
constant function and imposes exactness on every row.

For the dimensionless norm `Ĥ` and the skew part `A` of the weighted
derivative, solve `A*X = Ĥ*(h*Y) - B*X/2`. The allowed skew entries have
an optional bandwidth and are grouped under reflection. If `L` maps those
entries to `vec(A*X)`, the norm-only constraint is that the right-hand
side be orthogonal to `ker(L')`. This includes the symmetric compatibility
identity of [Glaubitz, Nordström and Öffner](https://arxiv.org/abs/2203.05479v2)
and also accounts for restrictions on the allowed skew entries.

After checking norm positivity, solve for the skew entries and check the
nullspace of the weighted derivative. Try the available free directions
when the particular solution has extra null modes. A failed candidate
search leaves feasibility unresolved.

The regression cases check exact algebraic identities for linear polynomials
on a translated interval with both norm structures. A numerical example fits
sine and cosine at all seven grid nodes using a bandwidth-two weighted
derivative. The tests also check reflection, rank, and rejection of
inconsistent value/derivative samples.
