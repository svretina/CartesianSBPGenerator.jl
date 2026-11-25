# Equal-order obstruction for local centered stencils

Let the repeating antisymmetric interior stencil have exact even order
`p`, and suppose the boundary closure is exact on all polynomials through
degree `p`. Work in dimensionless coordinates. Contract the symmetric
norm compatibility equation with the coefficients of the two-variable
polynomial `g(x,y) = (y-x)^p`.

The norm side vanishes because `(∂x + ∂y)g = 0`. The endpoint term also
vanishes because `g(0,0) = 0`. Since `p` is even, the two coupling
terms agree. The remaining compatibility equation is

```text
0 = 2 sum_{i,j across the cut} C[i,j] * (j-i)^p
  = 2 sum_{ℓ=1}^s c[ℓ] * ℓ^(p+1).
```

The second equality counts the `ℓ` pairs crossing the cut at offset
`ℓ`, using a closure at least as deep as the stencil radius `s`.
The final sum is the full interior moment of degree `p+1`, because
the stencil is antisymmetric. That moment is nonzero for a stencil of
exact order `p`. If it vanishes, the even moment at degree `p+2` is
also zero by antisymmetry, giving interior accuracy of at least `p+2`.

This contradiction occurs before checking norm positivity. It does
not depend on the closure depth or on the entries permitted in the
symmetric norm block. It applies to the local repeating-stencil ansatz.
The fixed-grid function-space method does not prescribe that interior
stencil and is a different construction.

The default order policy reports the obstruction immediately.
`order_policy=:allow_upgrade` explicitly permits increasing the
interior order while retaining the requested boundary order. The
returned family records both the original request and the achieved order.
