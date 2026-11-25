# Analytical alternating SBP operators

`Experimental.analytical_perturbation(p; N)` constructs a diagonal-norm
first derivative with two alternating interior stencils and exact order
$p$ throughout the grid. The checked orders are 2, 4, and 6. The method
starts with a conventional SBP(2p,p) family and uses exact rational
arithmetic. You can also pass an exact conventional family in place of
`p`, and supply a rational `epsilon` or a `domain` keyword.

## Construction

Let the conventional operator satisfy $Q_0=HD_0$ and
$Q_0+Q_0^T=B$. Take the degree-$p+1$ forward difference coefficients
$d_i=(-1)^{p+1-i}\binom{p+1}{i}$ for $i=0,\ldots,p+1$, and form
the vectors $u=(d,0)$ and $v=(0,d)$. The tile $A=uv^T-vu^T$ is
skew-symmetric. Both $u$ and $v$ are orthogonal to samples of every
polynomial through degree $p$ on the tile, so $A$ annihilates those
samples.

Embed the tile starting at each grid point where it fits, with sign
$(-1)^{j-1}$ at the one-based start index $j$. Their sum $K$ still
annihilates polynomials through degree $p$, and
$Q=Q_0+\epsilon K$ has the same symmetric part $B$. Consequently,
$D=H^{-1}Q$ preserves both the SBP identity and accuracy through
degree $p$. For even $p$ and odd $N$, reflection reverses the
sign of each tile and preserves the parity of its start index. This
gives $D=-JDJ$, with $J$ the reversal matrix.

In a fully overlapped interior row $i$, the increment to its
dimensionless tap $hD_{i,i+k}$ is
$\epsilon(-1)^{i-1}\sum_a(-1)^{a-1}A_{a,a+k}$, with the sum over
valid one-based tile indices. This gives the two phase stencils without
solving a new system for the closure. The perturbation has half-bandwidth
$p+2$. The complete derivative also includes the conventional closure.
For $p=2$, the phase A taps at offsets $-4,\ldots,4$ are

```text
 1/192   0   7/96   -2/3   0   2/3   -7/96   0   -1/192
```

Phase B reverses the perturbation's sign. The default amplitude
`1/(12*2^(2p))` gives signed leading bulk coefficients of magnitude
$1/12$, $1/48$, and $1/192$ at orders 2, 4, and 6, respectively.
The constructor checks the bulk taps, reflection, positive norm,
row-wise exact order, and a one-dimensional
nullspace on the requested grid. The last two checks apply to the
chosen amplitude and baseline.

## Matched-grid comparison

The benchmark compares the original direct alternating seeds, the
analytical perturbation, and its unperturbed conventional baseline on
$[0,1]$ using 41, 61, and 81 nodes. The reported derivative
coefficients are dimensionless. The SAT spectral radius is
$h \rho(L)$ for unit-speed inflow advection. The RK4 CFL divides the
smaller of the advection and reflecting wave spectral step limits by
$h$.

All operators had exact nullity one at all three resolutions.
The table gives the 81-node results.

| Order | Construction | Closure depth | Half-bandwidth | Max. absolute hD | Norm condition | SAT radius | RK4 CFL |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | Direct | 5 | 4 | 5.27 | 2.82 | 7.89 | 0.364 |
| 2 | Analytical | 4 | 4 | 1.75 | 3.47 | 1.37 | 2.06 |
| 2 | Conventional | 4 | 3 | 1.74 | 3.47 | 1.37 | 2.06 |
| 4 | Direct | 9 | 8 | 1.75e3 | 4.72 | 2.03e3 | 1.40e-3 |
| 4 | Analytical | 8 | 7 | 2.26 | 6.98 | 1.73 | 1.64 |
| 4 | Conventional | 8 | 7 | 2.26 | 6.98 | 1.73 | 1.64 |
| 6 | Direct | 15 | 14 | 9.76e7 | 19.0 | 5.12e7 | 5.52e-8 |
| 6 | Analytical | 14 | 13 | 60.9 | 820 | 12.8 | 0.221 |
| 6 | Conventional | 14 | 13 | 59.7 | 820 | 13.4 | 0.210 |

The perturbation reduces the largest coefficients and increases the
spectral step limit relative to these direct seeds at all three orders.
The sixth-order analytical operator inherits the conventional norm's
condition number and its 13-node half-bandwidth. The perturbation
itself reaches eight nodes on either side.

The table below reports maximum-norm semidiscrete errors using the
matrix exponential at time 8 for sinusoidal inflow advection and
time 1 for a reflecting wave pulse. All three entries per cell are
ordered by increasing resolution.

| Order | Construction | Advection error | Wave error |
| ---: | --- | --- | --- |
| 2 | Direct | 4.40e-3, 1.39e-3, 5.86e-4 | 9.80e-2, 5.42e-2, 3.95e-2 |
| 2 | Analytical | 1.79e-3, 7.38e-4, 3.99e-4 | 9.97e-2, 6.18e-2, 4.35e-2 |
| 2 | Conventional | 1.25e-3, 3.61e-4, 1.55e-4 | 1.02e-1, 6.17e-2, 4.37e-2 |
| 4 | Direct | 4.06e-4, 7.16e-5, 2.22e-5 | 6.70e-2, 2.91e-2, 2.01e-2 |
| 4 | Analytical | 9.33e-5, 1.23e-5, 2.92e-6 | 5.55e-2, 2.58e-2, 1.78e-2 |
| 4 | Conventional | 1.16e-4, 1.50e-5, 3.56e-6 | 5.54e-2, 2.53e-2, 1.79e-2 |
| 6 | Direct | 6.43e-5, 4.93e-6, 1.61e-6 | 1.03e-1, 5.53e-2, 3.65e-2 |
| 6 | Analytical | 6.07e-5, 4.47e-6, 6.16e-7 | 1.06e-1, 1.03e-1, 4.00e-2 |
| 6 | Conventional | 5.28e-5, 3.47e-6, 4.65e-7 | 2.28e-1, 8.04e-2, 4.35e-2 |

The sixth-order direct seed forces a very small common RK4 step on
these grids. Float64 errors in the powered RK4 update persist when
halving that step, so the tabulated errors use the semidiscrete
evolution. The script also prints errors for smooth test functions,
the signed leading Taylor moments, and the full PDE diagnostics for
every construction and resolution.

To reproduce the diagnostics and plots, run from the package root with
Julia and Python with NumPy and Matplotlib.

```sh
julia --project=. experiments/compare_perturbation.jl --pde --matrix-dir=/tmp/sbp-perturbation
python3 experiments/plot_perturbation.py /tmp/sbp-perturbation derivations/figures
```

Each figure shows the sparsity pattern, a shared signed logarithmic
coefficient scale, endpoint and interior stencils, the norm weights,
and the signed coefficient of $h^p f^{(p+1)}$ across the grid.

- [Order 2](figures/alternating_locality_p2.png)
- [Order 4](figures/alternating_locality_p4.png)
- [Order 6](figures/alternating_locality_p6.png)
