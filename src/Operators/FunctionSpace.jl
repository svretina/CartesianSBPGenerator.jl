function _reflection_orbits(entries, n)
    return [i + j == n + 1 ? [(i, j)] : [(i, j), (n + 1 - j, n + 1 - i)]
            for (i, j) in entries if i + j <= n + 1]
end

function _function_space_norm(::Type{T}, n, norm, closure_depth) where {T}
    fixed = zeros(T, n, n)
    if norm === :diagonal || norm === :auto
        isnothing(closure_depth) ||
            throw(ArgumentError("closure_depth applies to the restricted-full function-space norm"))
        entries = [(i, i) for i in 1:n]
    elseif norm === :restricted_full
        m = something(closure_depth, 2)
        1 <= m <= (n - 1) ÷ 2 ||
            throw(ArgumentError("restricted-full norm needs two closure blocks and an interior node"))
        for i in (m + 1):(n - m)
            fixed[i, i] = one(T)
        end
        entries = [(1, 1); [(i, j) for j in 2:m for i in 2:j]]
    else
        throw(ArgumentError("norm must be :diagonal, :restricted_full, or :auto"))
    end
    orbits = _reflection_orbits(entries, n)
    return fixed, orbits
end

function _function_space_skew_map(X::Matrix{T}, orbits) where {T}
    columns = zeros(T, length(X), length(orbits))
    for (column, orbit) in enumerate(orbits)
        action = reshape(view(columns, :, column), size(X))
        for (i, j) in orbit
            action[i, :] .+= X[j, :]
            action[j, :] .-= X[i, :]
        end
    end
    return columns
end

function _function_space_matrix(fixed::Matrix{T}, orbits, entries, sign) where {T}
    result = copy(fixed)
    for (orbit, value) in zip(orbits, entries), (i, j) in orbit
        result[i, j] = value
        result[j, i] = sign * value
    end
    return result
end

function _rank_of_weighted_derivative(Q::Matrix{T}) where {T}
    if T <: Rational
        return _homogeneous_solution(Q, :exact).rank
    end
    return rank(Q; atol=1e-12, rtol=max(size(Q)...) * eps(T))
end

function _function_space_Q(Bhalf, orbits, solution)
    n = size(Bhalf, 1)
    parameters = evaluate(solution)
    Q = _function_space_matrix(Bhalf, orbits, parameters, -1)
    rank_Q = _rank_of_weighted_derivative(Q)
    for direction in eachcol(solution.nullspace)
        rank_Q == n - 1 && break
        for step in (1, -1)
            candidate = parameters + step * direction
            trial = _function_space_matrix(Bhalf, orbits, candidate, -1)
            trial_rank = _rank_of_weighted_derivative(trial)
            if rank_Q < trial_rank <= n - 1
                Q, rank_Q, parameters = trial, trial_rank, candidate
                break
            end
        end
    end
    rank_Q == n - 1 || error("function-space skew search found no constant-only nullspace")
    return Q
end

"""
    generate_sbp(basis; method=:function_space, N, domain=(0, 1), arithmetic=:exact,
                 norm=:diagonal, closure_depth=nothing, bandwidth=N-1)

Construct on a fixed uniform grid. `BasisSamples` gives values and physical
derivatives at the grid nodes and determines `N` by its row count. For a
`PolynomialBasis` or `MonomialBasis`, supply `N` to sample the basis.
The bandwidth bounds the distance of nonzero entries of `Q` from its diagonal.
"""
function generate_sbp(
    basis::Union{BasisSamples,PolynomialBasis,MonomialBasis};
    method::Symbol=:function_space,
    N::Union{Nothing,Integer}=nothing,
    domain::Tuple{<:Real,<:Real}=(0, 1),
    derivative_order::Integer=1,
    arithmetic::Symbol=:exact,
    norm::Symbol=:diagonal,
    closure_depth::Union{Nothing,Integer}=nothing,
    bandwidth::Union{Nothing,Integer}=nothing,
)
    method === :function_space || throw(ArgumentError("sampled spaces use method=:function_space"))
    derivative_order == 1 || throw(ArgumentError("function-space construction supports derivative_order=1"))
    arithmetic in (:exact, :numerical) ||
        throw(ArgumentError("arithmetic must be :exact or :numerical"))
    n = isnothing(N) && basis isa BasisSamples ? size(basis.X, 1) : N
    !isnothing(n) && n >= 2 || throw(ArgumentError("specify at least two grid nodes"))
    width = something(bandwidth, n - 1)
    1 <= width < n || throw(ArgumentError("bandwidth must be between 1 and N-1"))
    T = arithmetic === :exact ? BigRational : Float64
    convert_sample = arithmetic === :exact ? _exact : Float64
    a, b = domain
    isfinite(a) && isfinite(b) && b > a || throw(ArgumentError("domain endpoints must be finite and ordered"))
    a, b = convert_sample(a), convert_sample(b)
    h = (b - a) / (n - 1)
    isfinite(h) && h > 0 || throw(ArgumentError("invalid grid spacing in selected arithmetic"))
    x = collect(range(a, b; length=n))
    all(>(0), diff(x)) || throw(ArgumentError("grid nodes must be distinct in selected arithmetic"))
    samples = basis isa BasisSamples ? basis : sample_basis(basis, x; arithmetic=T)
    size(samples.X, 1) == n || throw(DimensionMismatch("sample rows must match the number of grid nodes"))
    all(isfinite, samples.X) && all(isfinite, samples.Y) ||
        throw(ArgumentError("basis values and derivatives must be finite"))
    X = hcat(ones(T, n), Matrix{T}(convert_sample.(samples.X)))
    Y = hcat(zeros(T, n), Matrix{T}(convert_sample.(samples.Y)))

    B = zeros(T, n, n)
    B[1, 1], B[n, n] = -one(T), one(T)
    Bhalf = B / 2
    skew_entries = [(i, j) for j in 2:n for i in max(1, j - width):(j - 1)]
    skew_orbits = _reflection_orbits(skew_entries, n)
    L = _function_space_skew_map(X, skew_orbits)
    test_rows = transpose(_homogeneous_solution(Matrix(transpose(L)), arithmetic).nullspace)
    Hfixed, norm_orbits = _function_space_norm(T, n, norm, closure_depth)
    scaled_Y = h .* Y
    norm_columns = map(norm_orbits) do orbit
        action = _function_space_matrix(zeros(T, n, n), [orbit], T[1], 1) * scaled_Y
        test_rows * vec(action)
    end
    norm_system = LinearConstraintSystem(
        hcat(norm_columns...), test_rows * vec(Bhalf * X - Hfixed * scaled_Y),
    )
    norm_solution = solve_constraints(norm_system; backend=arithmetic)
    isnothing(norm_solution.particular) &&
        error("function-space norm compatibility returned $(norm_solution.status)")
    weights = evaluate(norm_solution)
    Hhat = _function_space_matrix(Hfixed, norm_orbits, weights, 1)
    if !_positive_definite(Hhat) && nullity(norm_solution) > 0
        target = T[i != j ? 0 : i == 1 ? 1//2 : 1 for (i, j) in first.(norm_orbits)]
        metric = T[length(orbit) * (first(orbit)[1] == first(orbit)[2] ? 1 : 2)
                   for orbit in norm_orbits]
        weights = _near_norm_solution(norm_solution, target, metric, arithmetic)
        Hhat = _function_space_matrix(Hfixed, norm_orbits, weights, 1)
    end
    _positive_definite(Hhat) || error("function-space norm selection found no positive-definite candidate")

    skew_solution = solve_constraints(
        LinearConstraintSystem(L, vec(Hhat * scaled_Y - Bhalf * X)); backend=arithmetic,
    )
    isnothing(skew_solution.particular) &&
        error("function-space skew reconstruction returned $(skew_solution.status)")
    Q = _function_space_Q(Bhalf, skew_orbits, skew_solution)
    H = h .* Hhat
    D = H \ Q
    _closure_residual_ok(Q + transpose(Q), B) && _closure_residual_ok(D * X, Y) ||
        error("function-space operator failed the SBP identity or basis exactness check")
    diagnostics = (; method, norm=norm === :auto ? :diagonal : norm, bandwidth=width,
                   norm_rank=norm_solution.rank, norm_nullity=nullity(norm_solution),
                   skew_nullity=nullity(skew_solution), nullity=1, exact=arithmetic === :exact)
    return FunctionSpaceSBPOperator(1, x, D, H, B, Q, copy(samples.labels), diagnostics)
end
