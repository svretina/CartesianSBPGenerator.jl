# **Implementation specification — CartesianSBPGenerator.jl**

## **1. Scope and scientific objective**

Develop a standalone Julia package for the automatic construction of **first-derivative Cartesian summation-by-parts finite-difference operators** through linear algebra.

The package must generate operators from prescribed interior and boundary accuracy orders, rather than retrieve predefined coefficients.

The implementation should primarily follow three papers:

- **Albin & Klarmann (2016):** reduced linear compatibility systems and algorithmic existence searches.
- **Glaubitz, Nordström & Öffner (2023):** general function-space SBP construction.
- **Glaubitz et al. (2026), Why Summation by Parts Is Not Enough:** nullspace consistency and additional operator-quality requirements.

The purpose is to develop reusable mathematical infrastructure for the subsequent construction of spherical and hyperboloidal SBP operators, particularly those requiring high-order accuracy at the origin and future null infinity.

### **Current scope**

| Capability                                             | Include           |
| ------------------------------------------------------ | ----------------- |
| First-derivative SBP operators                         | Yes               |
| Independently specified interior and boundary accuracy | Yes               |
| Uniform Cartesian grids                                | Yes               |
| Diagonal norms                                         | Yes               |
| Restricted full norms                                  | Yes               |
| Automatic closure-width determination                  | Yes               |
| Exact rational construction                            | Yes               |
| LinearSolve.jl backend                                 | Yes               |
| Reduced norm compatibility method                      | Yes               |
| General accuracy-basis abstraction                     | Yes               |
| Additional higher-degree accuracy constraints          | Yes               |
| Nullspace consistency                                  | Yes               |
| Spectral diagnostics                                   | Yes               |
| Optimization-based operator construction               | No                |
| Truncation-error optimization                          | No                |
| Nonuniform grids                                       | No                |
| Second and higher derivative SBP operators             | Future            |
| Spherical and hyperboloidal evolution                  | Separate packages |

The construction must remain based on linear systems. Positivity verification and, where necessary, positivity-feasibility searches are permitted, but the package should not implement the optimization-based SBP construction methodology at this stage.

# **2. Mandatory literature and public repositories**

The agent must inspect the following papers and their accompanying implementations before writing the construction algorithms.

## **2.1. Albin & Klarmann — Algorithmic exploration of SBP operators**

**Paper:** Nathan Albin and Joshua Klarmann, *An Algorithmic Exploration of the Existence of High-Order Summation by Parts Operators with Diagonal Norm* (2016).

- [ArXiv abstract](https://arxiv.org/abs/1403.5750?utm_source=chatgpt.com)
- [Full paper PDF](https://arxiv.org/pdf/1403.5750?utm_source=chatgpt.com)

This is the principal mathematical reference for the package.

The agent must study the block decomposition, elimination of skew-symmetric unknowns, reduced norm compatibility equations, positive-norm existence tests and automatic exploration of closure dimensions.

The paper parameterizes operators through interior order, boundary order and closure dimension, and demonstrates the construction of very high-order diagonal-norm operators.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://arxiv.org\&sz=32)

arXiv+1



I did not find a verified standalone public implementation belonging to the authors. There is, however, a relevant discussion in the Julia SBP ecosystem:

[SummationByPartsOperators.jl — Albin & Klarmann coefficients, Issue #130](https://github.com/ranocha/SummationByPartsOperators.jl/issues/130?utm_source=chatgpt.com).

The agent should inspect this discussion and its linked resources.

## **2.2. Glaubitz, Nordström & Öffner — Function-space SBP**

**Paper:** Jan Glaubitz, Jan Nordström and Philipp Öffner, *Summation-by-Parts Operators for General Function Spaces*, SIAM Journal on Numerical Analysis, 61(2), 733–754 (2023).

- [ArXiv abstract](https://arxiv.org/abs/2203.05479?utm_source=chatgpt.com)
- [Full paper PDF](https://arxiv.org/pdf/2203.05479?utm_source=chatgpt.com)
- [Published article](https://doi.org/10.1137/22M1470141?utm_source=chatgpt.com)

**Official accompanying repository:**

[**jglaubitz/FSBP**](https://github.com/jglaubitz/FSBP?utm_source=chatgpt.com)

Original reproducibility code for the function-space SBP paper.

Copy clone command

The repository is explicitly associated with the paper.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://github.com\&sz=32)

GitHub+1



The agent should inspect the implementation of:

- Function-space basis evaluation.
- Derivative evaluation of the basis.
- Positive quadrature construction.
- Skew-symmetric operator reconstruction.
- Accuracy and SBP verification.

Only the linear-algebra-based components relevant to our uniform-grid construction should be adopted.

Do not introduce general optimization or node-location optimization algorithms.

## **2.3. Glaubitz, Iske, Lampert & Öffner — Why SBP is not enough**

**Paper:** Jan Glaubitz, Armin Iske, Joshua Lampert and Philipp Öffner, *Why Summation by Parts Is Not Enough* (2026).

- [ArXiv abstract](https://arxiv.org/abs/2602.10786?utm_source=chatgpt.com)
- [Full paper PDF](https://arxiv.org/pdf/2602.10786?utm_source=chatgpt.com)

**Official reproducibility repository:**

[**JoshuaLampert/2026_SBP_not_enough**](https://github.com/JoshuaLampert/2026_SBP_not_enough?utm_source=chatgpt.com)

Reproducibility experiments and operator-quality investigations.

Copy clone command

The agent should use this repository primarily as a validation reference. Its purpose is to prevent us from accepting formally valid SBP operators that nevertheless have undesirable numerical properties.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://arxiv.org\&sz=32)

arXiv+1



In particular, investigate:

- Nullspace consistency.
- Spurious numerical null modes.
- SAT-modified eigenvalue properties.
- Accuracy beyond the prescribed function space.
- Operator conditioning.
- Numerical PDE validation.

Do not import its optimization-based operator-construction algorithms.

## **2.4. Additional Julia repositories**

These should also be inspected, primarily to reuse existing infrastructure and avoid unnecessary implementation work.

| Repository                                                                                                                     | Relevant components                                                             |
| ------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------- |
| [SummationByPartsOperatorsExtra.jl](https://github.com/JoshuaLampert/SummationByPartsOperatorsExtra.jl?utm_source=chatgpt.com) | Function-space SBP machinery, operator representations, validation              |
| [SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl?utm_source=chatgpt.com)                 | Established SBP coefficients, matrix interfaces, efficient operator application |
| [LinearSolve.jl](https://github.com/SciML/LinearSolve.jl?utm_source=chatgpt.com)                                               | Numerical solution of linear systems                                            |

`SummationByPartsOperatorsExtra.jl` is particularly relevant because it already implements functionality for function-space SBP operators in Julia.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://github.com\&sz=32)

GitHub



### **Instructions for inspecting external repositories**

Clone the reference repositories outside the new package's source tree:

```
mkdir -p references
cd references

git clone https://github.com/jglaubitz/FSBP.git
git clone https://github.com/JoshuaLampert/2026_SBP_not_enough.git
git clone https://github.com/JoshuaLampert/SummationByPartsOperatorsExtra.jl.git
git clone https://github.com/ranocha/SummationByPartsOperators.jl.git

cd ..
```

The agent must:

1. Read the papers and identify corresponding implementation files.
2. Produce a short report mapping the mathematical algorithms to the existing code.
3. Identify reusable routines, interfaces and tests.
4. Check each repository's software license before copying source code.
5. Preserve attribution and required license notices for adapted code.
6. Reimplement the mathematical algorithms independently where no compatible reuse license is available.
7. Record the inspected repository commit hashes.

Publicly accessible code is not automatically licensed for redistribution.

Use the external repositories as development references, not mandatory runtime dependencies, unless a dependency has a clear technical advantage.

Do not copy their complete optimization frameworks into this package.

# **3. Public API and operator specification**

Use the notation

\\[ \operatorname{SBP}^{(d)}(p_i,p_b), \\]

where

\\[ \begin{aligned} d&=\text{derivative order},\\\ p_i&=\text{interior accuracy order},\\\ p_b&=\text{boundary accuracy order}. \end{aligned} \\]

For the current release,

\\[ d=1. \\]

Keep `derivative_order` in the public API to facilitate future extensions, but explicitly reject unsupported values.

### **Minimal user interface**

```
using CartesianSBPGenerator

op = generate_sbp(
    4, 2;
    derivative_order = 1,
    norm = :diagonal,
    N = 101,
    domain = (0.0, 1.0)
)

D = derivative_matrix(op)
H = mass_matrix(op)
Q = weighted_derivative_matrix(op)
```

### **Generate a reusable operator family**

```
family = design_sbp(
    4, 4;
    derivative_order = 1,
    norm = :restricted_full
)

op = instantiate(
    family;
    N = 101,
    domain = (0.0, 1.0)
)
```

The second example is a construction request. Its feasibility must be determined by the package.

### **Supported norm types**

```
norm = :diagonal
norm = :restricted_full
norm = :auto
```

### **Advanced interface**

```
family = design_sbp(
    6, 6;
    derivative_order = 1,
    norm = :restricted_full,

    closure_depth = :auto,
    interior_radius = :minimal,

    construction = :reduced,
    arithmetic = :exact,

    extra_accuracy = :auto,

    max_closure_depth = 32,
    require_nullspace_consistency = true
)
```

No optimization-related keyword should be introduced in the initial implementation.

The requested accuracy orders must be respected. A request for `(4,4)` must not automatically return `(8,4)`.

# **4. Mathematical foundations**

## **4.1. Dimensionless formulation**

Consider a uniform boundary-inclusive grid,

\\[ x_i=a+ih, \qquad i=0,\ldots,N-1, \\]

with

\\[ h=\frac{b-a}{N-1}. \\]

Define dimensionless coordinates

\\[ \xi_i=i. \\]

Write the first-derivative operator and norm as

\\[ D=\frac1h\widehat D, \qquad H=h\widehat H. \\]

Then

\\[ HD=\widehat H\widehat D. \\]

The first-derivative SBP identity is

\\[ \boxed{ HD+D^TH=B, } \\]

where

\\[ B=e\_{N-1}e\_{N-1}^T-e_0e_0^T. \\]

Equivalently,

\\[ \widehat H\widehat D+ \widehat D^T\widehat H=B. \\]

Perform the entire coefficient construction in dimensionless coordinates.

Only introduce the physical grid spacing during instantiation.

## **4.2. Weighted derivative**

Define

\\[ Q=\widehat H\widehat D. \\]

The SBP identity reduces to

\\[ \boxed{Q+Q^T=B.} \\]

Parameterize

\\[ \boxed{ Q=\frac12B+A, \qquad A=-A^T. } \\]

Consequently,

\\[ \widehat D=\widehat H^{-1}Q. \\]

The central requirements of the package are therefore

\\[ \boxed{ \begin{aligned} Q+Q^T&=B,\\\ QX&=\widehat HY,\\\ \widehat H&=\widehat H^T\succ0. \end{aligned}} \\]

Here \\(X\\) and \\(Y\\) encode the function-space accuracy conditions.

The first two conditions are linear in the unknown coefficients.

Positive definiteness is checked separately.

# **5. Interior stencil construction**

For even interior accuracy

\\[ p_i=2s, \\]

use the standard minimum-width centered stencil,

\\[ (\widehat Df)\_i = \sum\_{j=-s}^{s}c_jf\_{i+j}. \\]

Enforce antisymmetry,

\\[ c\_{-j}=-c_j, \qquad c_0=0. \\]

The moment conditions are

\\[ \sum\_{j=-s}^{s}c_jj^k = \delta\_{k1}, \qquad k=0,\ldots,2s. \\]

By antisymmetry, only the odd moments need to be imposed.

Thus, solve

\\[ \boxed{ 2\sum\_{j=1}^{s}c_jj^{2\ell+1} = \delta\_{\ell0}, \qquad \ell=0,\ldots,s-1. } \tag{1} \\]

This is a square linear system of size \\(s\\).

### **Reference stencils**

\\[ p_i=2: \qquad \frac12[-1,0,1], \\]

\\[ p_i=4: \qquad \frac1{12}[1,-8,0,8,-1], \\]

\\[ p_i=6: \qquad \frac1{60}[-1,9,-45,0,45,-9,1]. \\]

The general minimum-width solution is

\\[ \boxed{ c_j= \frac{(-1)^{j+1}}{j} \frac{(s!)^2}{(s-j)!(s+j)!}. } \tag{2} \\]

Use Eq. (2) only as an independent validation formula. The production stencil generator must derive coefficients by solving Eq. (1).

### **Implementation requirements**

```
stencil = centered_stencil(
    interior_order = 6,
    radius = 3,
    arithmetic = Rational{BigInt}
)
```

Return:

```
stencil.offsets
stencil.coefficients
stencil.radius
stencil.accuracy_order
stencil.leading_error
```

Initially, use minimum-width interior stencils. Wider stencils may be supported as an explicitly requested linear construction experiment, but must not trigger dispersion or other optimization algorithms.

# **6. Boundary and norm structures**

Let

- \\(s\\): centered stencil radius.
- \\(m\\): number of closure rows at each boundary.
- \\(N\\): total number of grid points.

Require \\(m\geq s\\), with sufficient grid size that the left and right closures do not interact.

Partition the domain into

\\[ I_L=\\{0,\ldots,m-1\\}, \\]

\\[ I_C=\\{m,\ldots,N-m-1\\}, \\]

\\[ I_R=\\{N-m,\ldots,N-1\\}. \\]

## **6.1. Diagonal norm**

Parameterize the left norm block as

\\[ H_L= \operatorname{diag} (w_0,\ldots,w\_{m-1}). \\]

The full dimensionless norm is

\\[ \widehat H= \begin{pmatrix} H_L&0&0\\\ 0&I&0\\\ 0&0&R_mH_LR_m \end{pmatrix}. \\]

Here \\(R_m\\) is the reversal matrix.

There are \\(m\\) independent norm unknowns.

## **6.2. Restricted full norm**

Parameterize

\\[ \boxed{ H_L= \begin{pmatrix} w_0&0\\\ 0&\widetilde H_L \end{pmatrix}, } \\]

where

\\[ \widetilde H_L=\widetilde H_L^T. \\]

The endpoint remains decoupled from the remaining norm block.

There are

\\[ n_H=1+\frac{m(m-1)}2 \\]

independent unknowns.

**Diagonal norm**

**Restricted full norm**

Do not implement unrestricted full norms or interior implicit norms.

A bandwidth restriction within the restricted full block can be exposed as an optional parameter if necessary.

# **7. Albin–Klarmann reduced construction**

This is the principal algorithm of the package.

The original paper develops the compatibility equations in Section 2 and the existence algorithm in Section 3. It explicitly separates solving the linear equality constraints from determining whether a positive norm exists.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://arxiv.org\&sz=32)

arXiv



## **7.1. Construct the boundary-to-interior coupling**

Write

\\[ Q= \begin{pmatrix} A_L-\frac12e_0e_0^T&C&0\\\ -C^T&Q_C&\cdots\\\ 0&\cdots&Q_R \end{pmatrix}, \\]

where

\\[ A_L=-A_L^T. \\]

The interior stencil uniquely determines \\(Q_C\\) and the coupling matrix \\(C\\).

For the left boundary,

\\[ C\_{ij}=c\_{j-i}, \\]

whenever the column \\(j\\) is an interior node and the centered stencil reaches boundary row \\(i\\). All other entries vanish.

The skew-symmetry relation fixes the transpose block.

For example, for the second-order interior stencil and \\(m=3\\),

\\[ C= \begin{pmatrix} 0&0&\cdots\\\ 0&0&\cdots\\\ \frac12&0&\cdots \end{pmatrix}. \\]

Only the boundary row adjacent to the interior couples to the first interior node.

For fourth- and sixth-order interior stencils, the coupling extends over two and three rows, respectively.

The agent should generate this coupling automatically from stencil offsets.

## **7.2. Define the boundary accuracy basis**

For boundary accuracy \\(p_b\\), use

\\[ \mathcal F_b= \operatorname{span} \\{1,\xi,\xi^2,\ldots,\xi^{p_b}\\}. \\]

Define

\\[ X\_{ik}=\xi_i^k, \qquad i\in I_L, \\]

\\[ Y\_{ik}=k\xi_i^{k-1}. \\]

The derivative of the constant column is zero.

Let \\(X_C\\) contain the same basis functions evaluated at the interior nodes reached by \\(C\\).

The boundary accuracy conditions are

\\[ \left( A_L-\frac12e_0e_0^T \right)X+CX_C=H_LY. \\]

Therefore,

\\[ \boxed{ A_LX=F(H_L), } \tag{3} \\]

where

\\[ \boxed{ F(H_L) = H_LY+\frac12e_0e_0^TX-CX_C. } \tag{4} \\]

The challenge is to find a positive-definite \\(H_L\\) such that Eq. (3) admits a skew-symmetric solution.

## **7.3. Eliminate the skew-symmetric unknowns**

Since

\\[ A_L^T=-A_L, \\]

we have

\\[ X^TA_LX+(X^TA_LX)^T=0. \\]

Using Eq. (3),

\\[ \boxed{ X^TF+F^TX=0. } \tag{5} \\]

Provided \\(X\\) has full column rank, this condition is necessary and sufficient for the existence of a skew-symmetric \\(A_L\\) satisfying Eq. (3).

Substitute Eq. (4):

\\[ \boxed{ X^TH_LY+Y^TH_LX = X^TCX_C+X_C^TC^TX -X^Te_0e_0^TX. } \tag{6} \\]

Define the known matrix

\\[ \mathcal B\_{\mathrm{eff}} = X^TCX_C+X_C^TC^TX -X^Te_0e_0^TX. \\]

The final reduced compatibility equation is

\\[ \boxed{ X^TH_LY+Y^TH_LX =\mathcal B\_{\mathrm{eff}}. } \tag{7} \\]

This equation is linear in the entries of \\(H_L\\).

**We do not initially solve for the derivative stencil. We solve for the mass matrix only.**

This is the exact conceptual correspondence with our spherical SBP construction.

## **7.4. Diagonal-norm equations in scalar form**

For a diagonal norm,

\\[ H_L=\operatorname{diag}(w_0,\ldots,w\_{m-1}). \\]

The \\((j,k)\\) entry of the left-hand side of Eq. (7) is

\\[ \begin{aligned} &\sum\_{i=0}^{m-1} w_i\xi_i^j k\xi_i^{k-1}\\\ &\quad+ \sum\_{i=0}^{m-1} w_i j\xi_i^{j-1}\xi_i^k. \end{aligned} \\]

Hence,

\\[ \boxed{ (j+k) \sum\_{i=0}^{m-1} w_i\xi_i^{j+k-1} = (\mathcal B\_{\mathrm{eff}})\_{jk}, } \tag{8} \\]

for \\(j+k>0\\).

This has an important implementation consequence.

For diagonal norms, all equations with the same \\(j+k\\) have the same left-hand side. Their right-hand sides must therefore agree.

If they do not, the candidate operator is algebraically incompatible.

This provides a cheap consistency screening rule before a full linear solve.

For restricted full norms, retain the complete bilinear form in Eq. (7).

## **7.5. Assemble the reduced linear system**

Parameterize

\\[ H_L(z) = H_0+\sum\_{\alpha=1}^{n_H}z\_\alpha H\_\alpha. \\]

Define

\\[ K\_\alpha = X^TH\_\alpha Y+Y^TH\_\alpha X. \\]

Equation (7) becomes

\\[ \sum\_\alpha z\_\alpha K\_\alpha = \mathcal B\_{\mathrm{eff}} -X^TH_0Y-Y^TH_0X. \\]

Since each matrix is symmetric, vectorize only its upper triangular entries.

Define

\\[ \operatorname{vech}(M) \\]

as the vector of independent entries of a symmetric matrix.

Then

\\[ \boxed{ \mathsf A\_{:\alpha} = \operatorname{vech}(K\_\alpha), } \\]

and

\\[ \boxed{ \mathsf b = \operatorname{vech} \left( \mathcal B\_{\mathrm{eff}} -X^TH_0Y-Y^TH_0X \right). } \\]

Solve

\\[ \boxed{\mathsf Az=\mathsf b.} \tag{9} \\]

### **Pseudocode**

```
function assemble_reduced_system(X, Y, XC, C, H0, Hbasis)
    m = size(X, 1)

    e0 = zeros(eltype(X), m)
    e0[1] = one(eltype(X))

    E0 = e0 * transpose(e0)

    Beff =
        transpose(X) * C * XC +
        transpose(XC) * transpose(C) * X -
        transpose(X) * E0 * X

    rhs =
        Beff -
        transpose(X) * H0 * Y -
        transpose(Y) * H0 * X

    A = hcat([
        vech(
            transpose(X) * Hj * Y +
            transpose(Y) * Hj * X
        )
        for Hj in Hbasis
    ]...)

    return A, vech(rhs)
end
```

Implement this generically for diagonal and restricted full norms.

The diagonal-norm implementation may exploit Eq. (8) for improved efficiency.

## **7.6. Rank conditions**

The reduced theorem assumes that

\\[ \operatorname{rank}(X)=p_b+1. \\]

For a polynomial basis on distinct nodes, this is guaranteed when

\\[ m\geq p_b+1. \\]

Some classical constructions use smaller closure dimensions. For those, the full-column-rank assumption does not hold.

The agent must support these cases using the general condition

\\[ \boxed{ \ker X\subseteq\ker F } \\]

in addition to skew compatibility, or fall back to the direct linear-system construction.

This is important for reproducing minimal standard operators such as SBP(2,1) with a one-row closure.

Do not apply the full-column-rank theorem outside its assumptions.

## **7.7. Reconstruct the skew-symmetric boundary block**

After obtaining a positive-definite norm \\(H_L\\), reconstruct

\\[ A_LX=F(H_L). \\]

There are two ways to accomplish this.

### **Method A — Direct linear reconstruction**

Parameterize

\\[ A_L=\sum\_{\alpha}a\_\alpha E\_\alpha, \\]

where

\\[ E\_{ij}=e_ie_j^T-e_je_i^T, \qquad i\<j. \\]

Then

\\[ \sum\_\alpha a\_\alpha E\_\alpha X=F. \\]

Vectorize the equation and solve the resulting linear system.

This is the recommended implementation for the exact rational backend.

### **Method B — Analytical reconstruction**

For full-column-rank \\(X\\), define its pseudoinverse \\(X^+\\) and the projector

\\[ P=XX^+. \\]

A particular skew-symmetric solution is

\\[ \boxed{ A_L^{(0)} = FX^+-(FX^+)^T -(X^+)^TX^TFX^+. } \\]

The complete solution family is

\\[ \boxed{ A_L=A_L^{(0)}+(I-P)Z(I-P), \qquad Z=-Z^T. } \\]

In exact arithmetic, the pseudoinverse expression can be evaluated using exact linear solves.

In floating-point arithmetic, use QR or SVD factorizations rather than explicitly forming the inverse of \\(X^TX\\).

If

\\[ k=p_b+1 \\]

and \\(X\\) has rank \\(k\\), the dimension of the remaining skew freedom is

\\[ \boxed{ n\_{\mathrm{free},A} = \frac{(m-k)(m-k-1)}2. } \\]

This assumes an unrestricted dense skew block within the closure.

The agent should verify this nullity formula against the direct linear reconstruction.

# **8. Direct construction backend**

Although the reduced construction should be the default, also implement a simultaneous construction of \\(H_L\\) and \\(A_L\\).

This backend serves three purposes:

- Independent verification of the reduced algorithm.
- Construction of cases with rank-deficient basis matrices.
- Imposition of additional sparsity and accuracy constraints.

Parameterize

\\[ H_L=H_0+\sum\_\alpha z\_\alpha H\_\alpha, \\]

\\[ A_L=A_0+\sum\_\beta a\_\beta E\_\beta. \\]

Using

\\[ A_LX-H_LY = \frac12e_0e_0^TX-CX_C, \\]

assemble

\\[ \boxed{ \begin{pmatrix} \mathsf A_H&\mathsf A_A \end{pmatrix} \begin{pmatrix} z\\\a \end{pmatrix} =\mathsf b. } \\]

All equations are linear.

### **Required comparison**

For identical norm and skew sparsity assumptions, the direct and reduced methods must describe equivalent feasible operator families.

This should be tested automatically.

```
reduced = solve_constraints(
    spec;
    construction = :reduced
)

direct = solve_constraints(
    spec;
    construction = :direct
)
```

Do not require both methods to return identical particular solutions when the solution is nonunique. Instead, verify that each particular solution belongs to the other method's affine solution space.

# **9. Function-space abstraction**

This section implements the essential ideas of the 2023 Glaubitz paper.

## **9.1. General SBP exactness**

Let

\\[ \mathcal F=\operatorname{span} \\{f_1,\ldots,f_k\\}. \\]

Define

\\[ X\_{ij}=f_j(x_i), \\]

\\[ Y\_{ij}=f_j'(x_i). \\]

An \\(\mathcal F\\)-exact first derivative satisfies

\\[ DX=Y. \\]

With \\(Q=HD\\),

\\[ \boxed{QX=HY.} \\]

The SBP compatibility condition implies

\\[ \boxed{ X^THY+Y^THX=X^TBX. } \tag{10} \\]

For full-column-rank \\(X\\), a symmetric positive-definite norm satisfying this equation admits a compatible skew-symmetric weighted derivative.

This is the global analogue of our reduced boundary compatibility equation.

## **9.2. Connection to quadrature**

For a diagonal norm,

\\[ H=\operatorname{diag}(w_0,\ldots,w\_{N-1}), \\]

Eq. (10) becomes

\\[ \begin{aligned} \sum\_{i=0}^{N-1}w_i \big[ f_j(x_i)f_k'(x_i) +f_j'(x_i)f_k(x_i) \big] ={}&f_j(b)f_k(b)\\\ &-f_j(a)f_k(a). \end{aligned} \\]

Equivalently,

\\[ \boxed{ \sum\_{i=0}^{N-1} w_i(f_jf_k)'(x_i) = \int_a^b(f_jf_k)'(x)\\,dx. } \\]

This establishes the connection between function-space SBP operators and positive quadrature rules.

For

\\[ \mathcal F=\mathcal P_p, \\]

the required product-derivative space is

\\[ (\mathcal P_p\mathcal P_p)' =\mathcal P\_{2p-1}. \\]

Therefore, a diagonal norm must integrate polynomials through degree \\(2p-1\\) exactly if the derivative is globally exact on \\(\mathcal P_p\\).

For restricted full norms, retain the bilinear compatibility equation rather than replacing it with a scalar quadrature condition.

These observations follow from the general FSBP framework.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://epubs.siam.org\&sz=32)

SIAM+1



## **9.3. Software abstraction**

The construction engine must accept arbitrary basis evaluation matrices.

```
struct BasisSamples{T,MX,MY}
    X::MX
    Y::MY
    labels::Vector{String}
end
```

The concrete implementation should use appropriate parametric field types.

Provide

```
basis = PolynomialBasis(4)
```

and

```
basis = MonomialBasis([0, 1, 3, 5])
```

with

```
X = evaluate(basis, nodes)
Y = evaluate_derivative(basis, nodes)
```

Also provide an advanced interface:

```
basis = BasisSamples(
    X,
    Y,
    labels
)
```

This allows users to provide analytical basis evaluations directly.

### **Why this abstraction is essential for our project**

In the future spherical and hyperboloidal packages, we will need accuracy conditions on spaces such as

\\[ \mathcal F\_{\mathrm{scalar}} =\operatorname{span} \\{1,r^2,r^4,r^6\\}, \\]

\\[ \mathcal F\_{\mathrm{vector}} =\operatorname{span} \\{r,r^3,r^5\\}, \\]

and

\\[ \mathcal F\_{\mathscr I} =\operatorname{span} \\{x,x^2,x^3,\ldots\\}. \\]

We might also replace the ordinary derivative evaluations by

\\[ Y\_{ij} =(\mathcal L f_j)(x_i) \\]

for an analytically regularized differential operator \\(\mathcal L\\).

The current Cartesian package should not implement those operators, but its underlying constraint engine should accommodate such inputs without architectural changes.

Note that the spherical gradient–divergence pair requires a different compatibility identity from an ordinary first derivative. Reusing the basis and linear-system infrastructure does not imply that the Cartesian SBP theorem applies unchanged.

# **10. Linear solver requirements**

## **10.1. Numerical backend**

Use `LinearSolve.jl` for numerical linear systems.

```
using LinearAlgebra
using LinearSolve

problem = LinearProblem(A, b)

solution = solve(
    problem,
    QRFactorization(ColumnNorm())
)

z = solution.u
```

For poorly conditioned or rank-deficient systems:

```
solution = solve(
    problem,
    SVDFactorization()
)
```

These algorithms are documented by LinearSolve for systems where rank deficiency is possible.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://docs.sciml.ai\&sz=32)

LinearSolve.jl+1



The agent must verify the API against the pinned package version.

Do not assume that a successful solver return code implies an exact solution of the SBP compatibility equations.

Independently check the residual.

## **10.2. Exact backend**

For the polynomial uniform-grid construction, all coefficients are rational.

Use

```
const Exact = Rational{BigInt}
```

throughout the exact assembly pipeline.

Implement or adopt a suitable exact linear-algebra backend for:

- Row-reduced echelon form.
- Exact rank.
- Exact consistency.
- Particular solutions.
- Nullspace bases.

The agent should investigate existing Julia packages for exact rational linear algebra before implementing a custom elimination algorithm.

If no suitable lightweight dependency is available, implement exact Gaussian elimination with pivoting and rational arithmetic.

LinearSolve remains the numerical backend; do not assume its floating-point QR or SVD methods perform exact rational elimination.

## **10.3. Return an affine solution space**

For

\\[ \mathsf Az=\mathsf b, \\]

return

\\[ \boxed{ z=z_p+Z\theta, } \\]

where

\\[ \mathsf Az_p=\mathsf b, \\]

\\[ \mathsf AZ=0. \\]

Suggested output:

```
struct AffineSolution{T}
    particular::Vector{T}
    nullspace::Matrix{T}

    rank::Int
    augmented_rank::Int
    nullity::Int

    status::Symbol
end
```

The status must distinguish:

| Status                   | Interpretation                               |
| ------------------------ | -------------------------------------------- |
| `:unique`                | Exactly one solution                         |
| `:underdetermined`       | Consistent affine family                     |
| `:inconsistent`          | No solution                                  |
| `:numerically_uncertain` | Floating-point rank or residual is ambiguous |

Do not set arbitrary free parameters to zero and conclude that no positive norm exists if that particular solution is indefinite.

# **11. Norm positivity without operator optimization**

There is an unavoidable distinction between solving the linear equality constraints and enforcing positive definiteness.

Suppose

\\[ H_L(\theta) = H_p+\sum_j\theta_jH_j. \\]

The compatibility equations may admit infinitely many solutions, only some of which are positive definite.

We are not implementing optimization-based SBP construction, but we still need a feasibility procedure.

## **11.1. Diagonal norms**

The requirement is

\\[ w_i(\theta)>0. \\]

This defines a system of linear inequalities.

The original Albin–Klarmann algorithm uses linear programming to determine whether positive norm weights exist. This is an existence test rather than an optimization of the derivative operator.

Implement the following sequence:

1. Test the particular solution.
2. If positive, accept it.
3. Otherwise, solve the linear feasibility problem for the remaining parameters.
4. Independently certify the positivity and exactness of the selected solution.

Use a configurable strictly positive margin,

\\[ w_i\geq\epsilon. \\]

An optional LP backend is acceptable specifically for positivity feasibility.

## **11.2. Restricted full norms**

Require

\\[ w_0>0, \\]

\\[ \widetilde H_L\succ0. \\]

The particular solution should first be tested through Cholesky factorization.

If it is indefinite, positivity of the affine family remains unresolved.

For the initial implementation, support either:

- A deterministic search over rational free parameters.
- An optional semidefinite feasibility solver.

The latter would solve

\\[ H_L(\theta)\succeq\epsilon I. \\]

Do not implement spectral-radius or truncation-error optimization through SDP.

Importantly, failure of a finite parameter search does not prove that a positive-definite norm cannot exist.

Use the status

```
:positive_norm_not_found
```

rather than

```
:no_positive_norm_exists
```

unless an appropriate infeasibility certificate has been obtained.

## **11.3. Exact positive-definiteness certification**

For rational symmetric matrices, positive definiteness can be verified through an exact \\(LDL^T\\) factorization.

If

\\[ H_L=L\Lambda L^T, \\]

with nonsingular triangular \\(L\\), require

\\[ \Lambda\_{ii}>0 \\]

for every diagonal pivot.

This avoids relying solely on floating-point eigenvalues.

Also report the floating-point condition number,

\\[ \kappa_2(H_L), \\]

as a diagnostic. A positive but nearly singular norm may still be undesirable for numerical use.

# **12. Additional accuracy conditions**

This is a major design requirement inherited from your spherical SBP construction.

Once the baseline compatibility system has been solved, additional free parameters should be used to impose higher-degree accuracy conditions.

Do not introduce an objective function minimizing errors. Instead, introduce additional **linear equality constraints**.

## **12.1. Baseline conditions**

For boundary order \\(p_b\\), require

\\[ QX=HY \\]

for

\\[ X= \begin{pmatrix} \mathbf1& \boldsymbol\xi& \cdots& \boldsymbol\xi^{p_b} \end{pmatrix}. \\]

## **12.2. Next-degree condition**

The first additional monomial is

\\[ f\_{\mathrm{extra}}(\xi)=\xi^{p_b+1}. \\]

Its derivative is

\\[ f'\_{\mathrm{extra}}(\xi) =(p_b+1)\xi^{p_b}. \\]

The corresponding linear equality is

\\[ \boxed{ Q\boldsymbol\xi^{p_b+1} = (p_b+1)H\boldsymbol\xi^{p_b}. } \\]

For example:

| Boundary order | Mandatory degree | First additional degree |
| -------------- | ---------------- | ----------------------- |
| 1              | 1                | 2                       |
| 2              | 2                | 3                       |
| 3              | 3                | 4                       |
| 4              | 4                | 5                       |
| 5              | 5                | 6                       |
| 6              | 6                | 7                       |

## **12.3. Row-selective additional accuracy**

For diagonal norms, we can impose the next-degree condition at selected rows.

For boundary row \\(i\\),

\\[ (Q\mathbf f)\_i = (H\mathbf f')\_i. \\]

Because \\(H\\) is diagonal, this gives

\\[ (D\mathbf f)\_i=f'\_i. \\]

This is directly analogous to the additional accuracy conditions used in your spherical paper.

For restricted full norms, this implication is not generally valid at an individual row.

Instead,

\\[ D\mathbf f-\mathbf f' = H^{-1}(Q\mathbf f-H\mathbf f'). \\]

The inverse norm can mix weighted residuals between rows.

Therefore, distinguish two stages.

**Before fixing the norm:** impose additional exactness over the entire coupled norm block, or over isolated degrees of freedom such as the endpoint.

**After fixing the norm:** impose actual row-selective derivative accuracy using

\\[ e_i^TH^{-1}Q\mathbf f=f'\_i. \\]

Since \\(H\\) is now known, this is linear in the remaining skew-symmetric coefficients.

This is the preferred strategy for restricted full norms.

## **12.4. Rank-adaptive constraint addition**

Suppose the current system is

\\[ A_0z=b_0 \\]

with affine solution space

\\[ z=z_p+Z\theta. \\]

For a candidate additional equation,

\\[ a^Tz=\beta, \\]

substitute the affine representation:

\\[ \boxed{ a^TZ\theta=\beta-a^Tz_p. } \\]

Three cases occur:

| Condition                         | Action                                |
| --------------------------------- | ------------------------------------- |
| \\(a^TZ=0\\) and residual zero    | Constraint already satisfied          |
| \\(a^TZ=0\\) and residual nonzero | Reject as inconsistent                |
| \\(a^TZ\neq0\\)                   | Constraint reduces the affine freedom |

For the third case, solve the updated affine system and check that positive norms remain feasible.

### **Suggested deterministic policy**

Prioritize additional conditions in this order:

1. Lowest additional polynomial degree.
2. Endpoint itself.
3. Rows progressively farther from the endpoint.
4. Higher polynomial degrees if freedom remains.

Do not require all free parameters to be eliminated.

Some should be retained if necessary for positivity or nullspace consistency.

The complete sequence of accepted and rejected additional accuracy conditions must be included in the construction report.

# **13. Operator quality — Incorporating the 2026 paper**

The paper by Glaubitz, Iske, Lampert and Öffner is essential because it demonstrates that satisfying the formal SBP requirements does not necessarily produce a useful differentiation operator.

Their additional requirements concern locality, unresolved modes and accuracy outside the exactness space.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://arxiv.org\&sz=32)

arXiv+1



We will incorporate these as **construction constraints, acceptance criteria and diagnostics**, without adopting their optimization algorithms.

The accompanying repository now uses the revised paper title *Beyond the summation-by-parts property: nullspace consistency, sparsity, and regularization for FSBP operators*. The agent should use the latest revision of arXiv:2602.10786.

## **13.1. Nullspace consistency**

The first derivative must satisfy

\\[ D\mathbf1=0. \\]

For our nonperiodic Cartesian grid, we additionally require

\\[ \boxed{ \ker D=\operatorname{span}\\{\mathbf1\\}. } \\]

Equivalently,

\\[ \boxed{ \operatorname{rank}(D)=N-1. } \\]

Since \\(H\\) is invertible,

\\[ \ker D=\ker Q. \\]

Therefore, the kernel can be examined using the weighted derivative matrix directly.

This is advantageous for exact rational rank calculations.

### **Why it matters**

A formally consistent derivative may fail to differentiate additional discrete modes. Nullspace testing detects this class of defect.

The generator must not accept an operator with additional exact null modes as its default production result.

## **13.2. Repair through remaining skew freedom**

Suppose

\\[ A_L=A_L^{(0)}+Z_A(\theta), \\]

where every perturbation satisfies

\\[ Z_A^T=-Z_A, \\]

\\[ Z_AX=0. \\]

Such perturbations preserve the SBP identity and polynomial exactness.

If the first particular solution has an undesirable nullspace, investigate deterministic rational perturbations within this remaining affine space.

For example, test small integer or rational values of the free parameters and select the first candidate satisfying the mandatory structural conditions.

This is a feasibility search, not an optimization algorithm.

Do not perturb the derivative matrix directly in a way that violates the SBP property.

If no suitable candidate is found, widen the boundary closure or revisit the additional accuracy constraints.

## **13.3. Sparsity and locality**

The agent must distinguish three matrices:

\\[ H,\qquad Q,\qquad D=H^{-1}Q. \\]

For diagonal \\(H\\), the sparsity patterns of \\(Q\\) and \\(D\\) coincide.

For restricted full \\(H\\), the inversion of the local norm block can introduce additional nonzero entries within the boundary stencil.

Therefore, measure:

```
report.norm_bandwidth
report.weighted_derivative_bandwidth
report.derivative_bandwidth
report.boundary_stencil_width
report.nonzeros
```

The construction should maintain:

- A fixed repeating centered stencil in the bulk.
- Local boundary closures.
- No direct coupling between the two boundaries.
- No resolution-dependent global dense differentiation matrix.

For a fixed operator family, application cost must scale as

\\[ O(N). \\]

## **13.4. Accuracy beyond the mandatory basis**

Even though we will not minimize augmented-basis errors, we should measure them.

Define

\\[ X\_{\mathrm{aug}} = \begin{pmatrix} \boldsymbol\xi^{p_b+1}& \cdots& \boldsymbol\xi^{p_b+k} \end{pmatrix}. \\]

Calculate

\\[ R\_{\mathrm{aug}} = DX\_{\mathrm{aug}}-Y\_{\mathrm{aug}}. \\]

Report the residuals separately for:

- First boundary point.
- Remaining boundary closure.
- Interior.
- Transition rows.

For better numerical interpretation, use local Taylor moments.

At grid point \\(i\\),

\\[ \boxed{ E\_{i,k} = \frac1{k!} \sum_j\widehat D\_{ij}(j-i)^k -\delta\_{k1}. } \\]

Then

\\[ (Df)\_i-f'(x_i) = \sum\_{k=0}^{\infty} h^{k-1}E\_{i,k}f^{(k)}(x_i). \\]

These coefficients provide a direct characterization of the actual stencil accuracy and leading truncation error.

The agent must calculate and export them.

## **13.5. SAT-modified eigenvalue property**

The raw differentiation matrix does not constitute a complete PDE evolution operator.

For positive-speed advection, define

\\[ D\_{\mathrm{SAT}} = D+H^{-1}e_0e_0^T. \\]

The evolution equation is

\\[ u_t=-D\_{\mathrm{SAT}}u. \\]

The SBP identity implies

\\[ HD\_{\mathrm{SAT}} +D\_{\mathrm{SAT}}^TH = e_0e_0^T+e\_{N-1}e\_{N-1}^T. \\]

Thus the eigenvalues satisfy

\\[ \operatorname{Re}\lambda(D\_{\mathrm{SAT}})\geq0. \\]

A stronger property is

\\[ \boxed{ \operatorname{Re}\lambda(D\_{\mathrm{SAT}})>0. } \\]

The package should check this property independently.

Distinguish:

- Energy stability, which follows from the semidiscrete identity.
- Strict spectral decay, which requires additional conditions.
- Nullspace consistency of the unmodified derivative.

These properties are related but not equivalent.

Use the strict eigenvalue property as a reported quality criterion, with an option to require it during construction.

Do not implement spectral optimization.

# **14. Automatic closure search**

The generator should automatically determine a suitable boundary closure size.

Let

\\[ p_i=2s \\]

and let \\(p_b=t\\).

The Albin–Klarmann notation is

\\[ (s,t,m). \\]

Our public API uses

\\[ (p_i,p_b,d), \\]

while \\(m\\) is a construction parameter.

## **14.1. Search procedure**

**Input**

Interior order, boundary order, norm type.

Generate centered stencil and initialize closure depth.

Assemble reduced mass-matrix system.

Inconsistent

Increase closure depth.

Consistent

Check positive norm.

Reconstruct derivative and impose additional accuracy.

Check SBP, accuracy, locality and nullspace consistency.

**Return validated operator**

Otherwise continue the bounded search.

### **Pseudocode**

```
function design_sbp(p_int, p_bnd; config)

    validate_specification(p_int, p_bnd, config)

    stencil = generate_centered_stencil(p_int)

    for norm_type in requested_norms(config)

        for m in closure_candidates(stencil, config)

            spec = build_candidate(
                p_int,
                p_bnd,
                norm_type,
                stencil,
                m
            )

            system = assemble_reduced_system(spec)

            norm_family = solve_affine_system(system)

            if !isconsistent(norm_family)
                record_failure(spec, :inconsistent)
                continue
            end

            H = find_positive_norm(norm_family)

            if isnothing(H)
                record_failure(spec, :positivity_unresolved)
                continue
            end

            skew_family = reconstruct_skew_family(H, spec)

            enriched = add_accuracy_constraints(
                H,
                skew_family,
                spec
            )

            candidate = select_structurally_valid_member(
                H,
                enriched,
                spec
            )

            report = validate(candidate)

            if report.passed
                return candidate
            end

            record_failure(spec, report)

        end

    end

    return ConstructionFailure(...)
end
```

The code above specifies the algorithmic sequence, not production-ready Julia.

## **14.2. Feasibility screening**

For conventional diagonal-norm SBP first derivatives with the standard repeating interior stencil, apply the established restriction

\\[ p_i\geq2p_b. \\]

An exact-order request such as

```
design_sbp(4, 4; norm = :diagonal)
```

should fail with an explanatory diagnostic rather than initiate a long futile search.

For restricted full norms, do not assume the same restriction.

However, do not assume that increasing the closure size guarantees existence either.

## **14.3. Search results**

The package should distinguish:

```
:success
:theoretical_obstruction
:inconsistent_linear_system
:positive_norm_not_found
:certified_positivity_infeasible
:nullspace_inconsistent
:validation_failed
:search_limit_reached
```

An exhausted finite closure search must not be described as a proof of global nonexistence.

# **15. Worked example — SBP¹(2,1)**

This example should be the first end-to-end unit test.

Take the second-order centered stencil,

\\[ c\_{-1}=-\frac12, \qquad c_1=\frac12. \\]

Choose a one-point left closure,

\\[ m=1. \\]

The left norm block is

\\[ H_L=(w_0). \\]

The left weighted derivative row must have the structure

\\[ Q\_{L,:} = \begin{pmatrix} -\frac12&\frac12&0&\cdots \end{pmatrix}. \\]

The constant accuracy condition gives

\\[ Q\mathbf1=0, \\]

which is already satisfied.

The linear accuracy condition is

\\[ Q\boldsymbol\xi=H\mathbf1. \\]

At the first point,

\\[ -\frac12(0)+\frac12(1)=w_0. \\]

Therefore,

\\[ \boxed{w_0=\frac12.} \\]

The dimensionless first row becomes

\\[ \widehat D\_{0,:} = \frac1{w_0}Q\_{0,:}, \\]

so

\\[ \boxed{ \widehat D\_{0,:} = [-1,1,0,\ldots]. } \\]

The complete operator is

\\[ D=\frac1h \begin{pmatrix} -1&1&0&\cdots&0\\\ -\frac12&0&\frac12&\cdots&0\\\ 0&-\frac12&0&\ddots&0\\\ \vdots&\ddots&\ddots&\ddots&\frac12\\\ 0&\cdots&0&-1&1 \end{pmatrix}, \\]

with norm

\\[ H=h\operatorname{diag} \left( \frac12,1,\ldots,1,\frac12 \right). \\]

The generator must recover these coefficients without using hard-coded boundary weights.

Because \\(X\\) does not have full column rank for this minimal one-row closure, this test must exercise either the generalized reduced compatibility conditions or the direct backend.

# **16. Package architecture**

Keep the initial architecture relatively small.

```
CartesianSBPGenerator.jl/
├── Project.toml
├── README.md
├── LICENSE
├── AGENTS.md
├── ROADMAP.md
├── THIRD_PARTY_NOTICES.md
├── src/
│   ├── CartesianSBPGenerator.jl
│   ├── Types.jl
│   ├── Bases.jl
│   ├── Stencils.jl
│   ├── Norms.jl
│   ├── BoundaryBlocks.jl
│   ├── Reflection.jl
│   ├── Construction/
│   │   ├── Reduced.jl
│   │   ├── Direct.jl
│   │   ├── SkewReconstruction.jl
│   │   └── AccuracyEnrichment.jl
│   ├── Solvers/
│   │   ├── LinearSolveBackend.jl
│   │   ├── ExactBackend.jl
│   │   ├── AffineSolution.jl
│   │   └── Positivity.jl
│   ├── Search.jl
│   ├── Operators.jl
│   ├── Validation.jl
│   └── Serialization.jl
├── test/
│   ├── runtests.jl
│   ├── bases.jl
│   ├── stencils.jl
│   ├── reduced_construction.jl
│   ├── direct_construction.jl
│   ├── reconstruction.jl
│   ├── published_operators.jl
│   ├── positivity.jl
│   ├── enrichment.jl
│   ├── nullspace.jl
│   ├── convergence.jl
│   └── energy.jl
├── docs/
│   └── src/
├── derivations/
├── experiments/
├── benchmarks/
└── reports/
```

Avoid creating an elaborate framework for unsupported derivative orders.

Keep the `derivative_order` field in the specification, but only implement first-derivative construction.

## **16.1. Dependencies**

| Dependency                     | Purpose                               |
| ------------------------------ | ------------------------------------- |
| `LinearAlgebra`                | Matrix factorizations and diagnostics |
| `SparseArrays`                 | Sparse matrix construction            |
| `LinearSolve.jl`               | Numerical linear-system solution      |
| `Test`                         | Automated validation                  |
| `TOML`                         | Construction reports                  |
| `SummationByPartsOperators.jl` | Independent regression tests          |
| `BenchmarkTools.jl`            | Optional performance benchmarks       |
| `Documenter.jl`                | Documentation                         |

Optional feasibility dependencies should be implemented through package extensions where practical.

Do not add:

- General optimization libraries for stencil selection.
- Automatic differentiation.
- NonlinearSolve.jl.
- Symbolic differentiation frameworks.
- PDE solvers as mandatory package dependencies.

None is necessary for the initial construction algorithm.

# **17. Operator representation and matrix reconstruction**

Store the generated family using dimensionless coefficients.

```
struct SBPFamily{T}

    interior_order::Int
    boundary_order::Int
    derivative_order::Int

    interior_radius::Int
    closure_depth::Int

    norm_type::Symbol

    interior_coefficients::Vector{T}

    H_left::Matrix{T}
    Q_left::Matrix{T}
    coupling_left::Matrix{T}

end
```

The production representation should be concrete and type-stable.

For instantiation, construct

\\[ \widehat H= \operatorname{diag} (H_L,I,R_mH_LR_m). \\]

The right weighted derivative closure is obtained through antisymmetric reflection.

With \\(R\\) denoting reversal of the entire grid,

\\[ \boxed{ R\widehat HR=\widehat H, } \\]

\\[ \boxed{ RQR=-Q, } \\]

and consequently,

\\[ R\widehat DR=-\widehat D. \\]

The physical matrices are

\\[ H=h\widehat H, \\]

\\[ D=\frac1h\widehat H^{-1}Q. \\]

### **Required accessors**

```
grid(op)

derivative_matrix(op)

gradient_matrix(op)

mass_matrix(op)

weighted_derivative_matrix(op)

boundary_matrix(op)

left_closure(op)

interior_stencil(op)

construction_report(op)
```

Support

```
op * u

mul!(du, op, u)

Matrix(op)

sparse(op)
```

For restricted full norms, cache the factorization of the small norm block.

Do not construct a dense global inverse.

# **18. Validation requirements**

Validation must be performed independently of the equations used during construction.

## **18.1. Exact algebraic validation**

For rational coefficients, verify exactly:

\\[ Q+Q^T=B, \\]

\\[ QX=HY, \\]

\\[ H=H^T, \\]

\\[ D\mathbf1=0. \\]

Verify positive definiteness using exact factorization where possible.

## **18.2. Pointwise accuracy validation**

Compute the actual local Taylor moments

\\[ E\_{i,k}. \\]

Determine the achieved order of every row.

The output should include a vector such as

```
report.row_accuracy_orders
```

This is particularly useful when additional higher-degree conditions improve selected boundary points.

The interior and boundary orders should not be inferred merely from the requested construction constraints.

## **18.3. Convergence tests**

Use analytic functions such as

\\[ f(x)=e^x, \\]

\\[ f(x)=\sin(2\pi x), \\]

\\[ f(x)=\cos(x+x^2). \\]

For multiple resolutions, measure

\\[ \epsilon\_{\mathrm{int}} = \max\_{i\in I\_{\mathrm{int}}} |(Df)\_i-f'(x_i)|, \\]

and

\\[ \epsilon\_{\mathrm{bnd}} = \max\_{i\in I\_{\mathrm{bnd}}} |(Df)\_i-f'(x_i)|. \\]

Use at least four refinements where truncation error dominates roundoff.

Measure separately:

- Endpoint convergence.
- Remaining closure convergence.
- Interior convergence.
- Global maximum norm.
- Global weighted norm.

The global norm must not hide a reduced pointwise boundary order.

## **18.4. Energy identity test**

Consider

\\[ u_t=-cDu-cH^{-1}e_0e_0^Tu, \qquad c>0. \\]

The discrete energy is

\\[ E_h=\frac12u^THu. \\]

Using the SBP identity,

\\[ \boxed{ \dot E_h = -\frac c2 \left( u_0^2+u\_{N-1}^2 \right). } \\]

Test this identity on many arbitrary vectors \\(u\\).

A full time-dependent PDE solver is not required for this algebraic verification.

An optional experiment should evolve a smooth outgoing pulse using an appropriate Runge–Kutta method.

## **18.5. Nullspace tests**

For several grid sizes, verify

\\[ \operatorname{rank}(Q)=N-1. \\]

Include both odd and even \\(N\\).

Measure singular values to detect near-null modes in addition to exact kernel defects.

A useful report includes:

```
report.nullity
report.kernel_basis
report.smallest_nonzero_singular_value
```

## **18.6. Spectral diagnostics**

For the SAT-modified derivative, calculate:

\\[ \max\operatorname{Re}\lambda(-D\_{\mathrm{SAT}}), \\]

\\[ \min\operatorname{Re}\lambda(D\_{\mathrm{SAT}}), \\]

and

\\[ h\\,\rho(D\_{\mathrm{SAT}}). \\]

Check that the dimensionless spectral radius remains bounded under refinement for a fixed local family.

Do not optimize these quantities in the initial release.

# **19. Construction campaign**

The package should be validated on established operators before attempting novel equal-order families.

## **Campaign A — Diagonal norms**

| Interior order | Boundary order | Purpose                                    |
| -------------- | -------------- | ------------------------------------------ |
| 2              | 1              | Minimal reference                          |
| 4              | 2              | Classical reference                        |
| 6              | 3              | Classical reference                        |
| 8              | 4              | Higher-order reference                     |
| 10             | 5              | Higher-order existence                     |
| 12             | 6              | Practical sixth-order boundary alternative |

Reproduce selected closure-width and existence results from Albin and Klarmann.

In particular, compare the minimum closure size obtained through automatic search with their published results under matching assumptions.

## **Campaign B — Restricted full norms**

| Interior order | Boundary order | Purpose                                              |
| -------------- | -------------- | ---------------------------------------------------- |
| 2              | 2              | Equal-order feasibility                              |
| 4              | 3              | Classical full/block-norm reference where applicable |
| 4              | 4              | Equal-order feasibility                              |
| 6              | 5              | High-order reference where compatible                |
| 6              | 6              | Equal-order feasibility                              |

Only compare against published operators whose norm structure actually matches the restricted full ansatz.

Do not assume that an unrestricted full-norm operator satisfies our endpoint-decoupling requirement.

## **Campaign C — Automatic enrichment**

For every feasible family:

1. Construct the baseline operator.
2. Record the free norm and skew parameters.
3. Add the next-degree accuracy constraints.
4. Record accepted and rejected equations.
5. Reconstruct the derivative.
6. Check positivity and nullspace consistency.
7. Compare the leading truncation-error coefficients.

This campaign should demonstrate whether additional linear exactness conditions improve individual boundary points without sacrificing the original requirements.

### **Construction report template**

| Quantity                       | Required      |
| ------------------------------ | ------------- |
| Requested orders               | Yes           |
| Achieved orders                | Yes           |
| Norm structure                 | Yes           |
| Minimum closure found          | Yes           |
| Interior stencil radius        | Yes           |
| Norm-system rank               | Yes           |
| Norm-system nullity            | Yes           |
| Skew-system nullity            | Yes           |
| Additional conditions accepted | Yes           |
| Minimum norm eigenvalue        | Yes           |
| Condition number               | Yes           |
| Exact SBP residual             | Yes           |
| Polynomial residuals           | Yes           |
| Nullspace dimension            | Yes           |
| SAT spectral properties        | Yes           |
| Construction time              | Yes           |
| Failure reason                 | If applicable |

All results must be computed, not inferred from theoretical expectations.

# **20. Implementation milestones**

### **Development sequence**

1. **M01 — Literature and repository audit**

   Read the three essential papers, inspect public implementations and document reusable algorithms.
2. **M02 — Algebraic foundation**

   Implement basis evaluation, centered stencils, exact linear solves and affine solution representations.
3. **M03 — Albin–Klarmann reduced construction**

   Implement the norm-only compatibility equations and skew reconstruction.
4. **M04 — Norm families**

   Diagonal and restricted full norms, positivity and reference operator reproduction.
5. **M05 — Automatic closure search**

   Determine closure depth, classify failures and preserve requested accuracy.
6. **M06 — Additional linear accuracy**

   Rank-adaptive higher-degree exactness conditions.
7. **M07 — Operator-quality certification**

   Nullspace consistency, sparsity, actual accuracy and SAT spectra.
8. **M08 — Equal-order construction experiments**

   Investigate (2,2), (4,4) and (6,6), and validate feasible alternatives.
9. **M09 — Release**

   Public API, documentation, serialization and downstream integration.

## **M01 — Repository and literature audit**

**Deliverables:**

- `derivations/literature.md`
- `derivations/external_implementations.md`
- `THIRD_PARTY_NOTICES.md`

The agent should identify the mathematical equation corresponding to every reused algorithm.

In particular, I verified that the original FSBP repository contains the following MATLAB files:

| File                   | Investigation                    |
| ---------------------- | -------------------------------- |
| `compute_FSBP.m`       | Overall FSBP construction        |
| `compute_QF.m`         | Weighted derivative construction |
| `compute_LSQF.m`       | Quadrature-related construction  |
| `generate_span_poly.m` | Polynomial basis representation  |
| `generate_span.m`      | General function-space handling  |
| `solve_linear_SAT.m`   | PDE validation                   |

The repository uses the BSD-2-Clause license. The 2026 reproducibility repository and `SummationByPartsOperatorsExtra.jl` use MIT licenses.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://github.com\&sz=32)

GitHub+2



These licenses permit reuse subject to their attribution and notice requirements.

The agent should inspect these specific files rather than reverse-engineer the entire repositories unnecessarily.

## **M02 — Generic algebraic foundation**

**Acceptance criteria:**

- Generate standard centered stencils through at least order twelve.
- Verify exact Taylor moments.
- Support arbitrary sampled accuracy bases.
- Solve square, rectangular and rank-deficient linear systems.
- Return exact affine solution spaces.
- Correctly identify inconsistent systems.

## **M03 — Reduced construction**

**Acceptance criteria:**

- Implement Eq. (7).
- Implement both skew reconstruction methods.
- Support rank-deficient basis cases.
- Implement the direct reference construction.
- Demonstrate equivalence of direct and reduced approaches.

## **M04 — Norm families**

**Acceptance criteria:**

- Recover classical diagonal-norm operators.
- Implement restricted full norm parameterization.
- Verify exact positive definiteness.
- Reject invalid norm structures.
- Distinguish lack of a positive candidate from proven infeasibility.

## **M05 — Automatic closure search**

**Acceptance criteria:**

- Reproduce selected Albin–Klarmann existence results.
- Search without manually specifying closure depth.
- Preserve the requested accuracy orders.
- Return useful diagnostics when construction fails.

## **M06 — Additional accuracy**

**Acceptance criteria:**

- Add higher-degree conditions as linear equalities.
- Track rank changes.
- Preserve norm positivity.
- Verify actual derivative accuracy at selected rows.
- Do not implement truncation-error minimization.

## **M07 — Operator certification**

**Acceptance criteria:**

- Exact SBP identity.
- Exact mandatory polynomial conditions.
- Correct norm structure and positivity.
- Nullspace consistency.
- Local sparsity.
- Independent convergence tests.
- Semidiscrete energy balance.
- SAT eigenvalue diagnostics.

## **M08 — Research experiments**

Construct or investigate:

\\[ \operatorname{SBP}^{(1)}(2,2), \\]

\\[ \operatorname{SBP}^{(1)}(4,4), \\]

\\[ \operatorname{SBP}^{(1)}(6,6). \\]

Document mathematical obstructions or search limitations if any of these requests fail.

For comparison, provide

\\[ \operatorname{SBP}^{(1)}(4,2), \quad \operatorname{SBP}^{(1)}(8,4), \quad \operatorname{SBP}^{(1)}(12,6). \\]

The latter are useful for our downstream hyperboloidal project even if strictly equal-order local operators cannot be constructed.

## **M09 — Documentation and integration**

Provide a documented interface suitable for downstream use:

```
family = design_sbp(
    4, 4;
    derivative_order = 1,
    norm = :restricted_full
)

op = instantiate(family; N = 201)

G = gradient_matrix(op)
```

Export the exact dimensionless closure coefficients and norm blocks.

# **21. Standing instructions for the ChatGPT coding agent**

Place the following operational requirements in `AGENTS.md`.

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

# **22. Final mathematical specification**

The package should ultimately implement the following construction.

Given

\\[ (p_i,p_b,d=1), \\]

generate the centered interior stencil through

\\[ \sum_jc_jj^k=\delta\_{k1}. \\]

Construct the boundary polynomial or function-space matrices

\\[ X,\qquad Y,\qquad X_C. \\]

Solve the reduced norm compatibility system

\\[ \boxed{ X^TH_LY+Y^TH_LX = X^TCX_C+X_C^TC^TX-X^Te_0e_0^TX. } \\]

Find a positive-definite solution,

\\[ H_L\succ0. \\]

Reconstruct the skew-symmetric matrix,

\\[ \boxed{ A_LX = H_LY+\frac12e_0e_0^TX-CX_C. } \\]

Use available free parameters to impose additional linear accuracy conditions.

Assemble

\\[ Q+Q^T=B. \\]

Recover

\\[ \boxed{D=H^{-1}Q.} \\]

Finally, certify

\\[ \boxed{ \begin{aligned} H&\succ0,\\\ HD+D^TH&=B,\\\ DX&=Y,\\\ \ker D&=\operatorname{span}\\{\mathbf1\\}. \end{aligned}} \\]

This combines the three essential papers into a single coherent construction methodology:

| Paper                        | What we adopt                                                                           |
| ---------------------------- | --------------------------------------------------------------------------------------- |
| Albin–Klarmann               | Reduced mass-matrix linear systems and existence search                                 |
| Glaubitz–Nordström–Öffner    | General function-space accuracy abstraction                                             |
| Glaubitz–Iske–Lampert–Öffner | Structural validation beyond the formal SBP property                                    |
| Your spherical SBP paper     | Selective additional exactness conditions and analytical regularity-driven construction |

The immediate objective is not to produce an optimally tuned SBP operator. It is to produce **mathematically valid, reproducible, high-order Cartesian operators from linear systems**, with sufficient diagnostics that we can trust them when they become inputs to our spherical and hyperboloidal constructions.
