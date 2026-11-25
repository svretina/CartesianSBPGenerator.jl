> Historical roadmap. The [implementation specification](specification.md) supersedes this document and defines the current milestone sequence.

Agreed. We should generalize the package specification so that **derivative order, interior accuracy order and boundary accuracy order are three independent parameters**, and restrict the initial implementation to diagonal and restricted full norms.

This also means that the package should not be called a gradient generator internally, since it should eventually support higher derivatives. I would retain the name `CartesianSBPGenerator.jl`.

The following specification supersedes the corresponding sections of the previous roadmap. The remaining requirements concerning linear-system assembly, exact arithmetic, positivity, automatic closure search and validation remain applicable.

# **Revised specification — CartesianSBPGenerator.jl**

## **1. Operator notation and public API**

Define an operator family as

\\[ \boxed{\operatorname{SBP}^{(d)}(p_i,p_b)} \\]

where:

- \\(d\\): derivative order.
- \\(p_i\\): interior truncation accuracy order.
- \\(p_b\\): boundary truncation accuracy order.

Thus,

\\[ D^{(d)}f = \frac{d^df}{dx^d}+O(h^{p_i}) \\]

in the interior, and

\\[ D^{(d)}f = \frac{d^df}{dx^d}+O(h^{p_b}) \\]

in the boundary closure.

### **Public interface**

```
using CartesianSBPGenerator

op = generate_sbp(
    4, 2;
    derivative_order = 1,
    norm = :diagonal,
    N = 101,
    domain = (0.0, 1.0)
)
```

For a second derivative with the same accuracy orders:

```
op = generate_sbp(
    4, 2;
    derivative_order = 2,
    norm = :diagonal,
    N = 101
)
```

For a first derivative with uniform fourth-order accuracy:

```
op = generate_sbp(
    4, 4;
    derivative_order = 1,
    norm = :restricted_full,
    N = 101
)
```

The last example represents a construction request, not a guarantee that such an operator exists.

### **Interpretation of different requests**

| Request     | Derivative | Interior accuracy | Boundary accuracy |
| ----------- | ---------- | ----------------- | ----------------- |
| `SBP¹(2,1)` | First      | 2                 | 1                 |
| `SBP¹(4,4)` | First      | 4                 | 4                 |
| `SBP¹(6,5)` | First      | 6                 | 5                 |
| `SBP²(4,2)` | Second     | 4                 | 2                 |
| `SBP²(6,6)` | Second     | 6                 | 6                 |
| `SBP⁴(6,4)` | Fourth     | 6                 | 4                 |

I would make the positional arguments strictly represent accuracy orders, with the derivative order supplied as a keyword.

The package must also distinguish the requested accuracy orders from the achieved orders. A request for `(4,4)` must never silently produce `(8,4)` unless the user explicitly permits an interior-order upgrade.

## **2. Supported norm structures**

Restrict the package to two norm families.

### **2.1. Diagonal norm**

\\[ H=h\operatorname{diag} (w_0,\ldots,w\_{m-1},1,\ldots,1, w\_{m-1},\ldots,w_0). \\]

The independent unknowns are the \\(m\\) boundary weights.

Require

\\[ w_i>0. \\]

For conventional first-derivative diagonal-norm SBP operators, the necessary accuracy restriction is

\\[ \boxed{p_i\geq 2p_b.} \\]

This should be implemented as a theoretical screening rule before attempting the linear solve.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://epubs.siam.org\&sz=32)

SIAM+1

### **2.2. Restricted full norm**

Use

\\[ \boxed{ H=h \begin{pmatrix} H_L&0&0\\\ 0&I&0\\\ 0&0&R_mH_LR_m \end{pmatrix}, } \\]

where the left closure block has the structure

\\[ H_L= \begin{pmatrix} w_0&0\\\ 0&\widetilde H_L \end{pmatrix}, \\]

with

\\[ w_0>0,\qquad \widetilde H_L=\widetilde H_L^T\succ0. \\]

The boundary point itself remains decoupled from the remaining norm block. This is the defining restriction in Strand's restricted full-norm construction.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://www.sciencedirect.com\&sz=32)

ScienceDirect+1

**Diagonal norm**

Independent diagonal boundary weights.

**Restricted full norm**

Dense closure blocks, with isolated endpoint weights.

For a boundary closure of size \\(m\\), the number of independent norm unknowns is

\\[ \boxed{ n_H=1+\frac{m(m-1)}2. } \\]

The package should optionally support a bandwidth restriction within \\(\widetilde H_L\\), but this is a sparsity parameter of the restricted full norm rather than a separate norm family.

### **Supported norm arguments**

```
norm = :diagonal
norm = :restricted_full
norm = :auto
```

The `:auto` mode should search both norm families without silently modifying the requested accuracy orders.

Unrestricted full norms are outside the initial scope.

# **3. Generalized finite-difference stencil generator**

The interior stencil generator must now support arbitrary derivative orders.

Let

\\[ D^{(d)}=\frac1{h^d}\widehat D^{(d)}. \\]

For a centered stencil of radius \\(s\\),

\\[ (\widehat D^{(d)}f)\_i = \sum\_{j=-s}^{s}c_jf\_{i+j}. \\]

Taylor expansion gives

\\[ \boxed{ \sum\_{j=-s}^{s}c_jj^k = d!\\,\delta\_{kd}. } \\]

These equations must hold for

\\[ k=0,\ldots,d+p_i-1. \\]

This is the general accuracy condition for a derivative of order \\(d\\) with truncation accuracy \\(p_i\\).

## **3.1. Stencil symmetry**

The centered stencil must satisfy

\\[ \boxed{c\_{-j}=(-1)^d c_j.} \\]

Therefore:

- Odd derivatives use antisymmetric stencils.
- Even derivatives use symmetric stencils.

For odd derivatives,

\\[ c_0=0. \\]

For even derivatives, \\(c_0\\) is an independent coefficient.

For standard centered stencils on uniform grids, the requested interior accuracy should initially be even.

## **3.2. Minimum stencil width**

For even interior accuracy \\(p_i\\), the minimum centered stencil radius is

\\[ \boxed{ s\_{\min} = \left\lfloor\frac{d+1}{2}\right\rfloor +\frac{p_i}{2}-1. } \\]

The corresponding number of points is

\\[ n\_{\mathrm{stencil}}=2s\_{\min}+1. \\]

| Derivative order | Second-order accuracy | Fourth-order accuracy | Sixth-order accuracy |
| ---------------- | --------------------- | --------------------- | -------------------- |
| 1                | 3 points              | 5 points              | 7 points             |
| 2                | 3 points              | 5 points              | 7 points             |
| 3                | 5 points              | 7 points              | 9 points             |
| 4                | 5 points              | 7 points              | 9 points             |
| 5                | 7 points              | 9 points              | 11 points            |
| 6                | 7 points              | 9 points              | 11 points            |

These are minimum-width centered finite-difference stencils. The generator must also permit wider stencils with additional free parameters.

## **3.3. API**

```
stencil = centered_stencil(
    derivative_order = 2,
    accuracy_order = 6,
    radius = 3,
    arithmetic = Rational{BigInt}
)
```

For example, the sixth-order accurate second-derivative stencil is

\\[ \boxed{ \frac1{180h^2} [2,-27,270,-490,270,-27,2]. } \\]

The generator must recover these coefficients by solving the moment equations.

# **4. Generalized boundary accuracy conditions**

For a derivative of order \\(d\\),

\\[ \frac{d^d}{dx^d}x^k = \begin{cases} 0,&k\<d,\\\\[1mm] \dfrac{k!}{(k-d)!}x^{k-d},&k\geq d. \end{cases} \\]

Consequently, boundary accuracy of order \\(p_b\\) requires polynomial exactness through degree

\\[ \boxed{k\_{\max}=d+p_b-1.} \\]

This replaces the first-derivative-specific condition \\(k\_{\max}=p_b\\) used in the previous roadmap.

### **Required polynomial degrees**

| Derivative order | Boundary accuracy 2 | Boundary accuracy 4 | Boundary accuracy 6 |
| ---------------- | ------------------- | ------------------- | ------------------- |
| 1                | \\(0,\ldots,2\\)    | \\(0,\ldots,4\\)    | \\(0,\ldots,6\\)    |
| 2                | \\(0,\ldots,3\\)    | \\(0,\ldots,5\\)    | \\(0,\ldots,7\\)    |
| 3                | \\(0,\ldots,4\\)    | \\(0,\ldots,6\\)    | \\(0,\ldots,8\\)    |
| 4                | \\(0,\ldots,5\\)    | \\(0,\ldots,7\\)    | \\(0,\ldots,9\\)    |

Define the matrices

\\[ X= \begin{pmatrix} \mathbf1&\boldsymbol\xi& \boldsymbol\xi^2&\cdots& \boldsymbol\xi^{k\_{\max}} \end{pmatrix} \\]

and

\\[ X^{(d)} = \begin{pmatrix} \partial\_\xi^d\mathbf1& \partial\_\xi^d\boldsymbol\xi& \cdots& \partial\_\xi^d\boldsymbol\xi^{k\_{\max}} \end{pmatrix}. \\]

For the dimensionless weighted operator

\\[ Q_d=\widehat H\widehat D^{(d)}, \\]

the generic accuracy condition is

\\[ \boxed{Q_dX=\widehat HX^{(d)}.} \\]

This remains linear in the entries of \\(Q_d\\) and \\(\widehat H\\).

The same constraint-assembly infrastructure can therefore be used for every derivative order.

# **5. Different derivative orders require different SBP identities**

**This is the most important mathematical adjustment to the previous roadmap.**

We cannot simply take the first-derivative SBP construction and change the polynomial accuracy conditions to construct higher derivatives.

The discrete integration-by-parts identity itself depends on the derivative order.

## **5.1. First derivative**

The standard continuum identity is

\\[ \int u v_x\\,dx = [uv]-\int u_xv\\,dx. \\]

The discrete counterpart is

\\[ \boxed{ HD_1+D_1^TH=B. } \\]

Define

\\[ Q_1=HD_1=\frac12B+A, \\]

where

\\[ A=-A^T. \\]

The linear accuracy equations are

\\[ \boxed{Q_1X=HX'.} \\]

This remains exactly the construction in the original roadmap.

## **5.2. Second derivative**

The continuum identity is

\\[ \int u v\_{xx}\\,dx = [uv_x]-\int u_xv_x\\,dx. \\]

The corresponding standard SBP decomposition is

\\[ \boxed{ D_2=H^{-1}(-M+BS), } \\]

where:

- \\(H\\) is the positive-definite norm.
- \\(M=M^T\succeq0\\) is the stiffness matrix.
- \\(S\\) contains first-derivative approximations at the endpoints.

This is the construction used by Mattsson and Nordström for narrow-stencil second-derivative SBP operators.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://www.sciencedirect.com\&sz=32)

ScienceDirect+1

In dimensionless form,

\\[ H=h\widehat H, \\]

\\[ M=h^{-1}\widehat M, \\]

\\[ S=h^{-1}\widehat S. \\]

Consequently,

\\[ \widehat D_2 = \widehat H^{-1} (-\widehat M+B\widehat S). \\]

The accuracy conditions become

\\[ \boxed{ (-\widehat M+B\widehat S)X = \widehat HX''. } \\]

These equations remain linear if the unknowns are:

- Entries of \\(\widehat H\\).
- Independent symmetric entries of \\(\widehat M\\).
- Endpoint derivative coefficients of \\(\widehat S\\).

The positivity constraints are imposed after solving the affine system.

### **Accuracy conditions for the boundary derivative**

The first-derivative boundary approximation \\(S\\) has its own accuracy order, denoted \\(p_s\\).

Prescribe

\\[ (S\mathbf x^k)\_{\partial} = k\mathbf x^{k-1}\_{\partial}, \\]

for

\\[ k=0,\ldots,p_s. \\]

A standard initial choice is

\\[ p_s=p_b+1. \\]

The package should permit overriding this parameter.

### **Example: SBP²(4,2)**

The mandatory second-derivative conditions are

\\[ \begin{aligned} D_2\mathbf1&=0,\\\ D_2\mathbf x&=0,\\\ D_2\mathbf x^2&=2\mathbf1,\\\ D_2\mathbf x^3&=6\mathbf x. \end{aligned} \\]

With \\(p_s=3\\), the boundary derivative conditions are

\\[ \begin{aligned} S\mathbf1&=0,\\\ S\mathbf x&=\mathbf1,\\\ S\mathbf x^2&=2\mathbf x,\\\ S\mathbf x^3&=3\mathbf x^2, \end{aligned} \\]

evaluated at the two endpoints.

This defines a linear construction problem for the second derivative, rather than forming the square of an existing gradient matrix.

## **5.3. Compatible second derivatives**

For our eventual numerical relativity applications, I would also support construction of first- and second-derivative operators sharing the same norm.

Let \\(D_1\\) be a previously generated first derivative.

A compatible second derivative can be expressed as

\\[ \boxed{ D_2= H^{-1} \left( -D_1^THD_1-R+BS \right), } \\]

where

\\[ R=R^T\succeq0. \\]

If the endpoint derivative agrees with the endpoint rows of \\(D_1\\),

\\[ S\_{\partial}=(D_1)\_{\partial}, \\]

the operator is fully compatible.

These are established distinctions in the second-derivative SBP literature.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://doi.org\&sz=32)

DOI

Provide three construction modes:

```
compatibility = :independent
compatibility = :compatible
compatibility = :fully_compatible
```

For example:

```
D1 = generate_sbp(
    4, 2;
    derivative_order = 1,
    norm = :diagonal
)

D2 = generate_sbp(
    4, 2;
    derivative_order = 2,
    reference_operator = D1,
    compatibility = :compatible
)
```

When \\(D_1\\) and \\(H\\) are prescribed, the remaining construction equations for \\(R\\) and \\(S\\) are linear.

However, if we simultaneously optimize \\(D_1\\), \\(H\\) and \\(D_2\\), the compatibility condition is generally nonlinear. This should not be part of the initial linear solver.

## **5.4. Third and higher derivatives**

Repeated integration by parts produces

\\[ \boxed{ \begin{aligned} \int u\\,\partial_x^dv\\,dx ={}&(-1)^d\int v\\,\partial_x^du\\,dx\\\ &+ \left[ \sum\_{j=0}^{d-1} (-1)^j (\partial_x^ju) (\partial_x^{d-1-j}v) \right]\_{\partial\Omega}. \end{aligned}} \\]

This suggests a discrete identity of the form

\\[ HD_d-(-1)^dD_d^TH = \sum\_{j=0}^{d-1} (-1)^jS_j^TBS\_{d-1-j}, \\]

where \\(S_j\\) approximates the \\(j\\)-th derivative at the endpoints and \\(S_0=I\\).

For odd derivatives, the interior weighted operator has a skew-symmetric structure.

For even derivatives, an appropriate symmetric stiffness operator and sign-definiteness condition are additionally required.

Dedicated SBP operators for third and fourth derivatives have been developed by Mattsson.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://www.sciencedirect.com\&sz=32)

ScienceDirect+1

There is an important limitation: if all the endpoint operators \\(S_j\\) are unknown, their products can make the construction nonlinear.

Therefore, the initial implementation should support:

1. Arbitrary derivative orders in the finite-difference stencil generator.
2. First and second derivatives in the complete SBP construction engine.
3. Third and fourth derivatives through subsequent derivative-specific construction strategies.

The API should be generic from the beginning, but the package must explicitly reject derivative orders for which a validated SBP construction backend has not yet been implemented.

# **6. Revised mathematical architecture**

The construction should be divided into a derivative-independent core and derivative-specific strategies.

**User specification**

\\((p_i,p_b,d,\text{norm})\\)

**Generic stencil and polynomial generator**

**First derivative**

\\(Q+Q^T=B\\)

**Second derivative**

\\(Q=-M+BS\\)

**Higher derivative strategies**

Derivative-specific boundary identities

**Generic linear constraint engine**

\\(\mathsf Az=\mathsf b\\)

**LinearSolve.jl + exact arithmetic backend**

**Positivity, accuracy enrichment and validation**

**Validated SBP operator family**

The central constraint engine must not know which derivative is being constructed. It receives matrix parameterizations, polynomial test vectors and linear constraints.

# **7. Revised package types**

```
abstract type AbstractSBPStrategy end

struct FirstDerivativeSBP <: AbstractSBPStrategy end

struct SecondDerivativeSBP <: AbstractSBPStrategy end
```

Future extensions:

```
struct ThirdDerivativeSBP <: AbstractSBPStrategy end

struct FourthDerivativeSBP <: AbstractSBPStrategy end
```

Define the user specification:

```
Base.@kwdef struct SBPSpec
    interior_order::Int
    boundary_order::Int
    derivative_order::Int = 1

    norm::Symbol = :auto
    order_policy::Symbol = :exact

    interior_radius::Union{Nothing,Int} = nothing
    closure_depth::Union{Nothing,Int} = nothing
    norm_bandwidth::Union{Nothing,Int} = nothing

    boundary_derivative_order::Union{Nothing,Int} = nothing

    compatibility::Symbol = :independent
    extra_accuracy::Symbol = :auto

    symmetric::Bool = true
end
```

The exact type design can be improved by the agent using enums and parametric dispatch rather than relying excessively on symbols.

### **Operator output**

```
struct SBPOperator{T}
    derivative_order::Int

    requested_interior_order::Int
    requested_boundary_order::Int

    achieved_interior_order::Int
    achieved_boundary_order::Int

    x::Vector{T}

    D
    H
    B
    auxiliary_operators

    diagnostics
end
```

The `auxiliary_operators` field should contain derivative-specific data:

| Derivative        | Auxiliary matrices                        |
| ----------------- | ----------------------------------------- |
| First             | \\(Q\\)                                   |
| Second            | \\(M,S\\)                                 |
| Compatible second | \\(M,S,R,D_1\\)                           |
| Higher            | Appropriate boundary derivative operators |

The production implementation should use concrete parametric field types to ensure Julia type stability.

# **8. Revised public interface**

## **Generate and instantiate**

```
op = generate_sbp(
    6, 3;
    derivative_order = 1,
    norm = :diagonal,
    N = 201
)
```

## **Obtain derivative and norm matrices**

```
D = derivative_matrix(op)
H = mass_matrix(op)
B = boundary_matrix(op)
```

Retain convenience aliases:

```
gradient_matrix(op)  # derivative_order == 1
laplacian_matrix(op) # derivative_order == 2
```

The general accessor should always be `derivative_matrix`.

## **Generate a reusable operator family**

```
family = design_sbp(
    6, 5;
    derivative_order = 1,
    norm = :restricted_full
)

op1 = instantiate(family; N = 101)
op2 = instantiate(family; N = 201)
```

## **Investigate the linear system**

```
spec = SBPSpec(
    interior_order = 4,
    boundary_order = 4,
    derivative_order = 1,
    norm = :restricted_full
)

system = assemble_constraints(spec)

solution_family = solve_constraints(
    system;
    backend = :exact
)
```

## **Additional accuracy conditions**

```
family = design_sbp(
    6, 5;
    derivative_order = 1,
    norm = :restricted_full,
    extra_accuracy = :auto
)
```

For derivative order \\(d\\), the next polynomial degree after the mandatory boundary conditions is

\\[ \boxed{k\_{\mathrm{extra}}=d+p_b.} \\]

This is the generalization of the additional accuracy conditions from your spherical SBP construction.

# **9. Revised package architecture**

```
CartesianSBPGenerator.jl/
├── Project.toml
├── README.md
├── AGENTS.md
├── ROADMAP.md
├── src/
│   ├── CartesianSBPGenerator.jl
│   ├── Types.jl
│   ├── Stencils/
│   │   ├── Centered.jl
│   │   ├── Moments.jl
│   │   └── PolynomialBasis.jl
│   ├── Norms/
│   │   ├── Diagonal.jl
│   │   └── RestrictedFull.jl
│   ├── Strategies/
│   │   ├── AbstractStrategy.jl
│   │   ├── FirstDerivative.jl
│   │   ├── SecondDerivative.jl
│   │   └── Compatibility.jl
│   ├── Construction/
│   │   ├── Parameterization.jl
│   │   ├── BoundaryBlocks.jl
│   │   ├── AccuracyConstraints.jl
│   │   ├── SBPConstraints.jl
│   │   ├── Reflection.jl
│   │   └── Assembly.jl
│   ├── Solvers/
│   │   ├── LinearSolveBackend.jl
│   │   ├── ExactBackend.jl
│   │   ├── AffineFamily.jl
│   │   └── Rank.jl
│   ├── Optimization/
│   │   ├── Positivity.jl
│   │   ├── AccuracyEnrichment.jl
│   │   └── TruncationError.jl
│   ├── Search/
│   │   ├── Candidates.jl
│   │   └── AutomaticSearch.jl
│   ├── Operators/
│   │   ├── SBPFamily.jl
│   │   ├── Instantiation.jl
│   │   └── LinearAlgebraInterface.jl
│   └── Diagnostics/
│       ├── Accuracy.jl
│       ├── EnergyIdentities.jl
│       ├── Positivity.jl
│       └── Reports.jl
├── test/
│   ├── runtests.jl
│   ├── stencils.jl
│   ├── first_derivative.jl
│   ├── second_derivative.jl
│   ├── restricted_full.jl
│   ├── compatibility.jl
│   ├── equal_order.jl
│   └── convergence.jl
├── benchmarks/
├── experiments/
├── derivations/
└── docs/
```

Third- and fourth-derivative strategy files should only be added when their mathematical constructions are implemented.

# **10. Revised implementation milestones**

| Milestone | Objective                                          | Priority            |
| --------- | -------------------------------------------------- | ------------------- |
| M01       | Generic derivative stencil generator               | Essential           |
| M02       | Generic linear constraint engine                   | Essential           |
| M03       | Diagonal-norm first-derivative SBP                 | Essential           |
| M04       | Restricted-full-norm first-derivative SBP          | Essential           |
| M05       | Automatic closure search and accuracy enrichment   | Essential           |
| M06       | Equal-order first-derivative investigation         | Research            |
| M07       | Independent second-derivative SBP                  | Essential extension |
| M08       | Compatible and fully compatible second derivatives | Extension           |
| M09       | Third- and fourth-derivative SBP                   | Future extension    |
| M10       | Optimization, benchmarking and documentation       | Release             |

## **M01 — Generic finite-difference stencil construction**

Implement the moment equations for arbitrary derivative order.

Validate derivative orders \\(d=1,\ldots,6\\), with even accuracy orders up to twelve.

Do not restrict the stencil-generation infrastructure to the derivative orders currently supported by the SBP backends.

## **M02 — Generic algebraic engine**

Implement:

- Matrix parameterizations.
- Polynomial exactness constraints.
- Arbitrary linear constraints.
- Exact and numerical linear solves.
- Rank and nullspace extraction.
- Consistency classification.

This milestone must be independent of derivative order.

## **M03 — First derivative, diagonal norm**

Recover established families:

\\[ \operatorname{SBP}^{(1)}(2,1), \\]

\\[ \operatorname{SBP}^{(1)}(4,2), \\]

\\[ \operatorname{SBP}^{(1)}(6,3), \\]

\\[ \operatorname{SBP}^{(1)}(8,4). \\]

Verify that the generator correctly identifies the theoretical accuracy restriction.

## **M04 — First derivative, restricted full norm**

Reproduce published restricted full-norm families where available, particularly the fourth-order interior, third-order boundary construction of Strand.

Investigate higher-order restricted full norms without assuming that coefficients from unrestricted full-norm families remain valid under the endpoint-decoupling restriction.

## **M05 — Automated construction search**

Search closure depth, interior stencil radius and restricted norm bandwidth.

Implement additional higher-degree accuracy conditions and positive-definite norm selection.

The package should produce a useful operator using only:

```
generate_sbp(p_int, p_bnd; derivative_order = d)
```

when the requested family is supported and feasible within the search limits.

## **M06 — Equal-order first-derivative investigation**

This remains the immediate scientific objective for our hyperboloidal project.

| Requested family | Norm            | Alternative if infeasible |
| ---------------- | --------------- | ------------------------- |
| SBP¹(2,2)        | Restricted full | SBP¹(4,2)                 |
| SBP¹(4,4)        | Restricted full | SBP¹(8,4)                 |
| SBP¹(6,6)        | Restricted full | SBP¹(12,6)                |

Test existence using exact linear systems and positive-definiteness constraints.

Do not treat the alternative as the requested family; report it separately.

## **M07 — Independent second derivative**

Implement the narrow-stencil second-derivative construction

\\[ D_2=H^{-1}(-M+BS). \\]

Generate operators using both supported norm structures.

Require

\\[ H\succ0,\qquad M\succeq0. \\]

Validate against the Mattsson–Nordström constructions.

## **M08 — Compatible derivative families**

Implement first- and second-derivative operators sharing the same norm.

Verify the decomposition

\\[ M=D_1^THD_1+R, \qquad R\succeq0. \\]

Support both compatible and fully compatible families.

## **M09 — Higher derivatives**

Study and implement the dedicated third- and fourth-derivative SBP constructions from the literature.

This should be a separate milestone because the higher-derivative boundary operators and energy identities require additional theoretical work.

## **M10 — Final validation**

For each generated operator, independently verify:

| Property                | Requirement                         |
| ----------------------- | ----------------------------------- |
| Derivative order        | Correct                             |
| Interior accuracy       | At least requested                  |
| Boundary accuracy       | At least requested                  |
| Reflection symmetry     | Correct for derivative parity       |
| Norm structure          | Matches user specification          |
| Norm positivity         | Strictly positive definite          |
| Stiffness positivity    | Where applicable                    |
| SBP identity            | Exact in rational arithmetic        |
| Polynomial residuals    | Exact in rational arithmetic        |
| Operator application    | \\(O(N)\\) for fixed local closures |
| Resolution independence | Verified                            |
| Actual convergence      | Measured separately in all regions  |

# **11. Additional literature**

The previous references remain applicable. Add the following as mandatory reading for the higher-derivative milestones.

| Reference                                                                                                                                                                   | Application                                          |
| --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------- |
| [Strand (1994), SBP approximations for d/dx](https://doi.org/10.1006/jcph.1994.1005?utm_source=chatgpt.com)                                                                 | First derivative, diagonal and restricted full norms |
| [Mattsson and Nordström (2004), SBP approximations of second derivatives](https://doi.org/10.1016/j.jcp.2004.03.001?utm_source=chatgpt.com)                                 | Independent second derivatives                       |
| [Mattsson (2014), Diagonal-norm SBP operators for third and fourth derivatives](https://www.sciencedirect.com/science/article/pii/S0021999114004343?utm_source=chatgpt.com) | Higher derivative strategies                         |
| [Linders, Lundquist and Nordström (2018), On the Order of Accuracy of SBP Operators](https://doi.org/10.1137/17M1139333?utm_source=chatgpt.com)                             | Feasibility restrictions                             |
| [Svärd and Nordström (2014), Review of SBP schemes](https://arxiv.org/abs/1311.4984?utm_source=chatgpt.com)                                                                 | Mathematical foundations and stability               |

# **12. Final instructions to the implementation agent**

The package must treat the following as independent quantities:

\\[ \boxed{ \begin{aligned} d&=\text{derivative order},\\\ p_i&=\text{interior accuracy order},\\\ p_b&=\text{boundary accuracy order}. \end{aligned}} \\]

Only two norm families are required:

- Diagonal.
- Restricted full.

The mathematical construction must use derivative-specific SBP identities and the generic polynomial exactness condition

\\[ \boxed{ Q_dX=HX^{(d)}. } \\]

All linear equality constraints should be assembled into a common system

\\[ \mathsf Az=\mathsf b \\]

and solved using LinearSolve.jl for numerical exploration, with an exact rational backend for certification.

Do not construct higher derivatives merely by taking powers of the first-derivative matrix. Such compositions can be exposed as separate utilities, but they are not equivalent to independently constructed narrow-stencil SBP operators.

**For the immediate hyperboloidal project, complete milestones M01–M06 first.** These provide the Cartesian gradient generator we need for the subsequent spherical mass-matrix construction.

The derivative-order abstraction should be designed from the beginning so that adding second and higher derivatives later does not require restructuring the package. No actual equal-order coefficient systems have been solved in this specification, so their feasibility remains an explicit research target rather than an assumed capability.
