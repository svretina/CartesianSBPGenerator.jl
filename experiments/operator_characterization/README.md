# Cartesian SBP operator characterization

The pinned Julia 1.10 experiment builds the operators and compares their
uniform-grid refinements. It writes CSV tables, a Markdown report, and
CairoMakie figures. The underlying diagnostics also accept external matrices.

Run from the repository root.

~~~sh
bench=experiments/operator_characterization
julia +1.10 --project="$bench" -e 'using Pkg; Pkg.instantiate()'
julia +1.10 --project="$bench" "$bench/run_all.jl" "$bench/smoke.toml" "$bench/results-smoke"
~~~

The smoke configuration covers orders 2, 4, and 6, three analytical amplitudes,
two resolutions, and their even-grid counterparts. It exercises every PDE
model and the advanced spectral diagnostics. Run the full research campaign
with the following command.

~~~sh
julia +1.10 --project="$bench" "$bench/run_all.jl" "$bench/config.toml" "$bench/results-full"
~~~

The full configuration includes grids through 1,025 nodes and up to 100
periods or crossings. Dense SVD sampling and long RK4 evolutions are the most
expensive stages. Add `--no-pde` or `--no-plots` to select stages. Diagnostic
checkpoints are saved after every completed operator/grid pair. If a run is
interrupted, pass `--resume` with the same configuration, output directory,
Julia version, and PDE selection. Keep the output directory on disk, especially
when running across reboots.

~~~sh
julia +1.10 --project="$bench" "$bench/run_all.jl" "$bench/config.toml" "$bench/results-full" --resume
~~~

The launcher runs diagnostics and CairoMakie in separate processes. Rendering
uses one serial worker per operator, releasing its plotting memory on exit.
Regenerate the figures from completed results with the same experiment
environment. The renderer skips completed operators. Add `--force` to
re-render existing figures.

~~~sh
julia +1.10 --project="$bench" "$bench/generate_figures.jl" "$bench/results-smoke"
~~~

## Diagnostic conventions

- **Certification.** The gate checks the canonical endpoint matrix, the SBP
  residual, norm symmetry and positivity, reflection, and independently
  measured polynomial order. Exact rank is tested through the configured
  limit. Larger grids report the cutoff and the numerical rank when an SVD is
  scheduled.
- **Structure and accuracy.** Data include signed local Taylor coefficients,
  left-closure ABTE and MBTE, alternating-phase cancellation, dimensionless
  matrices, bandwidths, norm spectra and bilinear quadrature residuals.
  Exact operators use BigFloat for smooth-function differentiation.
  Convergence rates stay within each grid parity and stencil phase.
- **Spectra and dispersion.** The raw derivative and inflow SAT spectra are
  evaluated independently. The energy norm uses a Cholesky similarity
  transform. The energy-CFL interval is sampled. Bulk spectra come from the
  repeating norm and weighted derivative symbols, including compact schemes.
  Eigenvector overlap tracks Bloch branches. A closed-form test checks the
  analytical branches. PPW uses the connected low-frequency interval below
  the error threshold. Sampled sign changes identify candidate zeros.
- **PDEs and energy.** Models include exact-inflow and periodic advection,
  outgoing packets, absorbing and reflecting waves, and manufactured
  variable-speed advection in both split forms. Power uses the evolved
  boundary driver. The integrated balance uses trapezoidal quadrature over
  saved times. Packet diagnostics include post-exit residual energy and an
  estimate of incoming branch energies from windowed spectra. The windowed
  estimates are sensitive to spectral leakage. Rigorous frequency-domain
  scattering requires a separate calculation.
- **Time integration and cost.** RK4 uses the spectral stability interval
  and an order-dependent accuracy cap. The final field is checked against
  a matrix-exponential reference and the step is reduced when needed.
  Semidiscrete propagation is available when the RK4 step demand exceeds
  the cap. Only time-resolved explicit runs enter work–precision plots.
  The PDE timer covers the sparse assembled RHS evolution loop. The
  derivative application benchmark measures dense, sparse, row-stencil,
  and cached norm-factor application separately.
- **Coverage.** STD, RF, ALT, ANA and the existing compact orders 4 and 6
  are compared. ALT and ANA have odd finite-domain grids. ALT-B was
  removed from the construction scope. The status table records these cases
  and the configured compact resolution cap. The equal-order comparison
  uses minimum pointwise order. The cost comparison matches same-grid
  operators using measured application time or nonzero count.

## Data and figures

| Output | Contents |
| --- | --- |
| report.md and order_*.csv | Matched-order tables and coverage |
| summary.csv | Measured orders, coefficient hashes, Taylor errors, norm and spectral measures |
| status.csv | Unsupported grids, resource limits, and observed deficiencies |
| differentiation.csv | Smooth-function errors, regions, and refinement rates |
| evolution.csv | PDE errors, time checks, runtime, phase, and branch estimates |
| performance.csv and construction.csv | Application timings and construction stages |
| comparable_cost.csv | Nearest application cost and nonzero count |
| metadata.toml and config.toml | CPU, Julia and BLAS settings, manifest hash, and effective configuration |
| operators/ | Per-grid arrays, exact sparse entries, and individual run traces |
| figures/ | Operator panels and grouped comparisons |

Plots use a consistent family palette, line styles and markers. The default
width is 510 pt for a two-column layout. PDFs retain vector text and lines,
embed fonts, and rasterize dense layers at 360 dpi. The full configuration
exports 300 dpi PNG previews. Exact zeros on logarithmic axes sit at the
plotting floor. [CairoMakie export settings](https://docs.makie.org/stable/explanations/backends/cairomakie)

The optional desktop viewer reads the saved evolution data from an
environment containing GLMakie.

~~~sh
julia --project=/path/to/glmakie-environment -i "$bench/interactive_evolution.jl" "$bench/results-smoke"
~~~

## External matrices

Supply physical matrices and a uniform grid. The diagnostic module infers
the grid spacing and uses it to compute dimensionless quantities.

~~~julia
using CartesianSBPGenerator.Diagnostics

op = OperatorData(D, H, B, x; depth=4, period=1)
report = diagnose_operator(op; minimum_order=2)
report.algebra.certified
~~~
