struct FirstDerivativeSystem{T}
    constraints::LinearConstraintSystem{T,Matrix{T},Vector{T}}
    spec::SBPSpec
    closure_depth::Int
    closure_width::Int
    canonical_size::Int
    interior_stencil::CenteredStencil{T}
    norm::Symbol
    norm_entries::Vector{Tuple{Int,Int}}
    construction::Symbol
    X::Matrix{T}
    Y::Matrix{T}
    C::Matrix{T}
    K::Matrix{T}
end

_default_closure_depth(boundary_order::Int) = boundary_order == 1 ? 1 : 2 * boundary_order

function _closure_method(method, construction)
    isnothing(method) && return something(construction, :reduced)
    method in (:classical, :reduced) ||
        throw(ArgumentError("stencil families support method=:classical or :reduced; " *
                            "use generate_sbp(basis; method=:function_space) for a sampled space"))
    selected = method === :classical ? :direct : :reduced
    isnothing(construction) || construction === selected ||
        throw(ArgumentError("method and construction select different backends"))
    return selected
end

function _closure_configuration(spec::SBPSpec)
    _validate_spec(spec)
    strategy = _strategy(spec)
    strategy isa FirstDerivativeSBP || throw(ArgumentError(
        "second-derivative constraint assembly is not implemented yet",
    ))
    spec.compatibility === :independent || throw(ArgumentError(
        "first-derivative construction only supports compatibility=:independent",
    ))
    spec.symmetric || throw(ArgumentError("first-derivative construction requires symmetric=true"))
    diagonal_orders = spec.boundary_order <= spec.interior_order ÷ 2
    norm = if spec.norm === :auto
        diagonal_orders && isnothing(spec.norm_bandwidth) ? :diagonal : :restricted_full
    else
        spec.norm
    end
    (norm !== :diagonal || diagonal_orders) || throw(ArgumentError(
        "a first-derivative diagonal norm requires interior_order >= " *
        "2*boundary_order; requested ($(spec.interior_order), $(spec.boundary_order))",
    ))
    if !isnothing(spec.norm_bandwidth)
        norm === :restricted_full ||
            throw(ArgumentError("norm_bandwidth applies only to the restricted full norm"))
        spec.norm_bandwidth >= 0 || throw(ArgumentError("norm_bandwidth must be nonnegative"))
    end

    radius = isnothing(spec.interior_radius) ?
        minimum_centered_radius(1, spec.interior_order) : spec.interior_radius
    stencil = centered_stencil(
        derivative_order=1, accuracy_order=spec.interior_order, radius=radius,
    )
    nullity(stencil) == 0 || throw(ArgumentError(
        "first-derivative SBP construction currently requires a minimum-radius interior stencil",
    ))
    default_depth = norm === :diagonal ? _default_closure_depth(spec.boundary_order) :
                    spec.boundary_order + 2
    closure_depth = something(spec.closure_depth, max(radius, default_depth))
    closure_depth >= radius || throw(ArgumentError("closure_depth must be at least the interior radius"))
    closure_width = closure_depth + radius
    n = 2 * closure_width + 1
    norm_entries = [(i, i) for i in 1:closure_depth]
    if norm === :restricted_full
        bandwidth = something(spec.norm_bandwidth, closure_depth - 1)
        append!(norm_entries, [
            (i, j) for j in 3:closure_depth for i in max(2, j - bandwidth):(j - 1)
        ])
    end
    return stencil, norm, norm_entries, closure_depth, closure_width, n
end

function _closure_samples(spec, basis, width, ::Type{T}) where {T}
    nodes = 0:(width - 1)
    samples = sample_basis(PolynomialBasis(spec.boundary_order), nodes; arithmetic=T)
    isnothing(basis) && return samples.X, samples.Y
    extra = basis isa BasisSamples ? basis : sample_basis(basis, nodes; arithmetic=T)
    size(extra.X, 1) == width ||
        throw(DimensionMismatch("basis samples must cover the $width closure stencil nodes"))
    all(isfinite, extra.X) && all(isfinite, extra.Y) ||
        throw(ArgumentError("basis samples must be finite"))
    convert_sample = T === BigRational ? _exact : Float64
    return hcat(samples.X, convert_sample.(extra.X)), hcat(samples.Y, convert_sample.(extra.Y))
end

function _boundary_coupling(stencil::CenteredStencil{T}, m) where {T}
    radius = stencil.radius
    C = zeros(T, m, radius)
    c = coefficients(stencil)
    for j in 1:radius, i in (m + j - radius):m
        C[i, j] = c[m + j - i + radius + 1]
    end
    return C
end

function _norm_actions(Y::Matrix{T}, entries) where {T}
    return map(entries) do (i, j)
        action = zeros(T, size(Y))
        action[i, :] .= Y[j, :]
        if i != j
            action[j, :] .= Y[i, :]
        end
        action
    end
end

_skew_entries(m) = [(i, j) for j in 2:m for i in 1:(j - 1)]

function _skew_columns(X::Matrix{T}) where {T}
    m = size(X, 1)
    entries = _skew_entries(m)
    columns = zeros(T, length(X), length(entries))
    for (k, (i, j)) in enumerate(entries)
        action = reshape(view(columns, :, k), size(X))
        action[i, :] .= X[j, :]
        action[j, :] .= -X[i, :]
    end
    return columns
end

function _homogeneous_solution(A::Matrix{T}, backend) where {T}
    solution = solve_constraints(LinearConstraintSystem(A, zeros(T, size(A, 1))); backend)
    isnothing(solution.particular) &&
        error("could not resolve the basis rank ($(solution.status))")
    return solution
end

function _compatibility_vector(X, F, kernel)
    bilinear = transpose(X) * F
    symmetry = [bilinear[i, j] + bilinear[j, i]
                for j in axes(bilinear, 2) for i in 1:j]
    return isempty(kernel) ? symmetry : vcat(symmetry, vec(F * kernel))
end

"""
    assemble_constraints(spec; method=:reduced, basis=nothing, arithmetic=:exact)

Assemble norm-only compatibility equations or a simultaneous norm and skew
system (`method=:classical`, also selected by `construction=:direct`).
Optional basis functions augment the required
polynomials. Supply `BasisSamples` at the dimensionless nodes of the left
closure stencil when using an analytical basis.
"""
function assemble_constraints(
    spec::SBPSpec;
    construction::Union{Nothing,Symbol}=nothing,
    method::Union{Nothing,Symbol}=nothing,
    basis::Union{Nothing,BasisSamples,PolynomialBasis,MonomialBasis}=nothing,
    arithmetic::Symbol=:exact,
)
    construction = _closure_method(method, construction)
    construction in (:reduced, :direct) ||
        throw(ArgumentError("construction must be :reduced or :direct"))
    arithmetic in (:exact, :numerical) ||
        throw(ArgumentError("arithmetic must be :exact or :numerical"))
    stencil, norm, entries, m, width, n = _closure_configuration(spec)
    T = arithmetic === :exact ? BigRational : Float64
    typed_stencil = CenteredStencil(
        stencil.derivative_order, stencil.accuracy_order, stencil.radius, stencil.offsets,
        T.(stencil.particular), T.(stencil.nullspace),
    )
    values, derivatives = _closure_samples(spec, basis, width, T)
    X, Y = values[1:m, :], derivatives[1:m, :]
    C = _boundary_coupling(typed_stencil, m)
    K = -(C * values[(m + 1):end, :])
    K[1, :] .+= X[1, :] ./ 2
    actions = _norm_actions(Y, entries)
    names = [i == j ? Symbol(:w_, i - 1) : Symbol(:h_, i - 1, :_, j - 1)
             for (i, j) in entries]
    if construction === :reduced
        kernel = _homogeneous_solution(X, arithmetic).nullspace
        matrix = hcat([_compatibility_vector(X, action, kernel) for action in actions]...)
        rhs = -_compatibility_vector(X, K, kernel)
    else
        matrix = hcat([-vec(action) for action in actions]..., _skew_columns(X))
        rhs = vec(K)
        append!(names, [Symbol(:a_, i - 1, :_, j - 1) for (i, j) in _skew_entries(m)])
    end
    constraints = LinearConstraintSystem(matrix, rhs, names)
    return FirstDerivativeSystem(constraints, spec, m, width, n, typed_stencil,
                                 norm, entries, construction, X, Y, C, K)
end

_backend(::FirstDerivativeSystem{<:Rational}) = :exact
_backend(::FirstDerivativeSystem{<:AbstractFloat}) = :numerical

function solve_constraints(system::FirstDerivativeSystem; backend::Symbol=_backend(system), kwargs...)
    return solve_constraints(system.constraints; backend, kwargs...)
end

function _skew_block(entries, m, ::Type{T}) where {T}
    A = zeros(T, m, m)
    for (value, (i, j)) in zip(entries, _skew_entries(m))
        A[i, j], A[j, i] = value, -value
    end
    return A
end

function _analytical_skew(X::Matrix{T}, F::Matrix{T}, backend) where {T}
    sample_space = _homogeneous_solution(X, backend)
    entries = _skew_entries(size(X, 1))
    independent = sample_space.pivot_columns
    Xr, Fr = X[:, independent], F[:, independent]
    A = zeros(T, size(X, 1), size(X, 1))
    if !isempty(independent)
        Xt = transpose(Xr)
        Xplus = if T <: Rational
            (Xt * Xr) \ Matrix(Xt)
        else
            qr(Xr) \ Matrix{T}(I, size(X, 1), size(X, 1))
        end
        P = Fr * Xplus
        correction = transpose(Xplus) * Xt * Fr * Xplus
        A .= P - transpose(P) - correction
        A .= (A - transpose(A)) / 2
    end
    if !_closure_residual_ok(A * X, F)
        return AffineSolution(:inconsistent, 0, Int[], nothing, zeros(T, length(entries), 0))
    end
    Z = _homogeneous_solution(Matrix(transpose(X)), backend).nullspace
    directions = _skew_entries(size(Z, 2))
    basis = T[Z[i, k] * Z[j, l] - Z[i, l] * Z[j, k]
              for (i, j) in entries, (k, l) in directions]
    status = isempty(directions) ? :unique : :underdetermined
    return AffineSolution(status, length(entries) - length(directions), Int[],
                          T[A[i, j] for (i, j) in entries], basis)
end

_closure_residual_ok(actual::Matrix{<:Rational}, target) = actual == target
_closure_residual_ok(actual::Matrix{<:AbstractFloat}, target) =
    isapprox(actual, target; rtol=1e-10, atol=1e-12)

"""
    reconstruct_skew(system, H; method=:linear)

Return the affine space of independent upper-triangular skew entries for a
fixed closure norm `H`. Use `method=:analytical` for the pseudoinverse formula.
"""
function reconstruct_skew(system::FirstDerivativeSystem{T}, H::AbstractMatrix;
                          method::Symbol=:linear) where {T}
    size(H) == (system.closure_depth, system.closure_depth) ||
        throw(DimensionMismatch("norm must match the closure dimensions"))
    X, F = system.X, Matrix{T}(H) * system.Y + system.K
    method === :analytical && return _analytical_skew(X, F, _backend(system))
    method === :linear || throw(ArgumentError("method must be :linear or :analytical"))
    return solve_constraints(LinearConstraintSystem(_skew_columns(X), vec(F));
                             backend=_backend(system))
end
