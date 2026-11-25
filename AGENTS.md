# CartesianSBPGenerator.jl — Agent instructions

## Mission

Implement a standalone Julia package that constructs first-derivative Cartesian SBP operators by solving linear systems. Do not retrieve hard-coded operators as a substitute for construction.

## Scope

- Uniform, boundary-inclusive Cartesian grids.
- Independently specified interior and boundary accuracy orders.
- First derivatives only.
- Diagonal and restricted full norms.
- Exact rational and LinearSolve.jl numerical backends.
- Reduced Albin–Klarmann compatibility construction.
- Generic function-space accuracy constraints.
- Additional higher-degree linear accuracy conditions.
- Automatic closure search.
- Structural and numerical certification.

Do not implement optimization-based SBP construction, truncation-error minimization, spectral optimization, nonuniform grids, higher derivatives, spherical operators or hyperboloidal PDEs.

## Mandatory references

1. Albin and Klarmann, arXiv:1403.5750.
2. Glaubitz, Nordström and Öffner, arXiv:2203.05479.
3. Glaubitz, Iske, Lampert and Öffner, arXiv:2602.10786.

Inspect the accompanying public repositories before implementing the mathematical algorithms.

Reuse compatible code where appropriate, preserve license notices, cite the original work and record the upstream commit hashes.

## Mathematical construction

Use Q = HD and parameterize Q through its prescribed symmetric boundary contribution and skew-symmetric remainder.

Construct the reduced norm compatibility equations before solving for the skew matrix.

Preserve the affine solution space of every underdetermined linear system.

Check positive definiteness separately from linear consistency.

Reconstruct D only after obtaining a valid norm and weighted derivative.

Never explicitly invert the global mass matrix during operator application.

## Accuracy requirements

The positional input parameters represent interior and boundary accuracy, respectively. The derivative order is an independent keyword and must currently equal one.

Do not silently lower boundary accuracy or increase interior order.

Verify actual row-wise accuracy after reconstructing the derivative matrix.

Use additional polynomial exactness constraints to exploit available free parameters. Do not introduce error-minimization objectives.

For restricted full norms, account for the mixing introduced by the inverse norm when imposing row-specific accuracy.

## Validation requirements

Every production operator must satisfy:

- Exact SBP compatibility for exact arithmetic.
- Mandatory polynomial exactness.
- Positive-definite norm.
- Correct norm structure.
- Reflection symmetry.
- Local boundary closure.
- Nullspace consistency.
- Independently measured convergence.

Report SAT-modified spectral properties, norm conditioning and leading truncation-error coefficients.

A successful linear solve alone does not constitute a successful operator construction.

## Development workflow

Work through the milestones sequentially.

Before implementing a milestone:

1. Read the relevant mathematical derivation.
2. Inspect existing reusable implementations.
3. Write or update the mathematical specification.
4. Implement the smallest complete solution.
5. Add independent regression tests.
6. Execute the complete test suite.
7. Record results and remaining issues.

Do not move to a new milestone with unresolved regressions.

Do not invent numerical results.

If a mathematical system is inconsistent, report the obstruction or the limits of the explored ansatz. Do not weaken mandatory constraints without explicitly documenting and requesting approval for the change.

## Engineering requirements

Prefer simple, type-stable Julia implementations.

Keep construction separate from operator instantiation.

Use sparse matrices or local stencil representations for large grids.

Use exact rational arithmetic for coefficient certification.

Avoid unnecessary dependencies and premature optimization.

Maintain reproducible construction reports and exact coefficient serialization.

The package must remain standalone and usable by subsequent spherical and hyperboloidal SBP packages.
