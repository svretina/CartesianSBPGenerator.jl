const R = BigRational

struct StencilExperiment
    constraints::LinearConstraintSystem{R,Matrix{R},Vector{R}}
    order::Int
    period::Int
    closure_depth::Int
    stencil_radius::Int
    norm_bandwidth::Int
    closure_norm_bandwidth::Int
    norm::Symbol
    Hfixed::Matrix{R}
    B::Matrix{R}
    norm_orbits::Vector{Vector{Tuple{Int,Int}}}
    skew_orbits::Vector{Vector{Tuple{Int,Int}}}
end

function _boundary_matrix(n)
    B = zeros(R, n, n)
    B[1, 1], B[n, n] = -1, 1
    return B
end

# A connection that touches an interior row is determined by its distance
# and the phase of its upper-triangular row. Reflection identifies phases.
function _orbit_key(i, j, n, m, period)
    ir, jr = n + 1 - j, n + 1 - i
    if m < i <= n - m || m < j <= n - m
        return (1, j - i, min(mod(i - 1, period), mod(ir - 1, period)))
    end
    u, v = min((i, j), (ir, jr))
    # Number right closures by their reflected indices, independent of N.
    return (0, u, v)
end

function _orbits(entries, n, m, period)
    groups = Dict{NTuple{3,Int},Vector{Tuple{Int,Int}}}()
    indices = NTuple{3,Int}[]
    for (i, j) in entries
        key = _orbit_key(i, j, n, m, period)
        if !haskey(groups, key)
            groups[key] = Tuple{Int,Int}[]
            push!(indices, key)
        end
        push!(groups[key], (i, j))
    end
    return [groups[key] for key in indices]
end

function _actions(X, orbits, sign)
    L = zeros(R, length(X), length(orbits))
    for (k, orbit) in enumerate(orbits)
        action = reshape(view(L, :, k), size(X))
        for (i, j) in orbit
            action[i, :] .+= X[j, :]
            if i != j
                action[j, :] .+= sign .* X[i, :]
            end
        end
    end
    return L
end

"""
    assemble_stencil_experiment(order; N, period=1, closure_depth, stencil_radius,
                                norm=:periodic, norm_bandwidth=0,
                                closure_norm_bandwidth=norm_bandwidth)

Assemble exact, fixed-grid SBP equations through degree `order`.
`:periodic` uses a symmetric banded norm with periodic interior entries, with
bandwidth zero giving diagonal weights. `:identity_interior` fixes the
interior norm to identity and permits dense symmetric closure blocks.
"""
function assemble_stencil_experiment(
    order::Integer;
    N::Integer,
    period::Integer=1,
    closure_depth::Integer,
    stencil_radius::Integer,
    norm::Symbol=:periodic,
    norm_bandwidth::Integer=0,
    closure_norm_bandwidth::Integer=norm_bandwidth,
)
    order >= 1 || throw(ArgumentError("order must be positive"))
    period >= 1 || throw(ArgumentError("period must be positive"))
    1 <= stencil_radius <= closure_depth ||
        throw(ArgumentError("require 1 ≤ stencil_radius ≤ closure_depth"))
    0 <= norm_bandwidth <= closure_depth || throw(ArgumentError("invalid norm bandwidth"))
    0 <= closure_norm_bandwidth <= closure_depth ||
        throw(ArgumentError("invalid closure norm bandwidth"))
    N >= max(order + 1, 2 * closure_depth + period) ||
        throw(ArgumentError("N must fit the basis and at least one full interior period"))
    norm in (:identity_interior, :periodic) || throw(ArgumentError("unsupported experimental norm"))
    norm === :identity_interior && (norm_bandwidth != 0 || closure_norm_bandwidth != 0) &&
        throw(ArgumentError("norm bandwidths apply to norm=:periodic"))
    n, m = Int(N), Int(closure_depth)
    same_closure(i, j) = j <= m || i > n - m
    q_entries = [(i, j) for j in 2:n for i in 1:(j - 1)
                 if same_closure(i, j) || j - i <= stencil_radius]
    Hfixed = zeros(R, n, n)
    if norm === :identity_interior
        h_entries = [(i, j) for j in 1:n for i in 1:j if same_closure(i, j)]
        for i in (m + 1):(n - m)
            Hfixed[i, i] = 1
        end
    else
        h_entries = [(i, j) for j in 1:n for i in 1:j
                     if j - i <= (same_closure(i, j) ? closure_norm_bandwidth : norm_bandwidth)]
    end
    q_orbits = _orbits(q_entries, n, m, period)
    h_orbits = _orbits(h_entries, n, m, period)
    B = _boundary_matrix(n)
    X, Y = polynomial_matrices(0:(n - 1), order, 1)
    A = hcat(-_actions(Y, h_orbits, 1), _actions(X, q_orbits, -1))
    b = vec(Hfixed * Y - (B / 2) * X)
    if norm === :periodic && period == 1
        # The repeating compact norm stencil is normalized to unit row sum.
        sums = _actions(ones(R, n, 1), h_orbits, 1)
        A = vcat(A, hcat(sums[(m + 1):(m + 1), :], zeros(R, 1, length(q_orbits))))
        push!(b, one(R))
    end
    return StencilExperiment(LinearConstraintSystem(A, b), order, period, m, stencil_radius,
                             norm_bandwidth, closure_norm_bandwidth, norm, Hfixed, B,
                             h_orbits, q_orbits)
end

function _row_orders(D, x)
    orders = fill(-1, length(x))
    values = ones(R, length(x))
    derivatives = zeros(R, length(x))
    for degree in 0:length(x)
        residual = D * values - derivatives
        for i in eachindex(orders)
            if orders[i] == -1 && !iszero(residual[i])
                orders[i] = degree - 1
            end
        end
        all(>=(0), orders) && break
        derivatives .= (degree + 1) .* values
        values .*= x
    end
    return orders
end

function _operator(H, D, B, Q, x, order, method; metadata=(;))
    row_orders = _row_orders(D, x)
    exactness = minimum(row_orders) >= order
    exactness && Q == H * D && Q + transpose(Q) == B ||
        error("experimental operator failed exact validation")
    count_nonzero(A) = count(!iszero, A)
    bandwidth(A) = maximum(abs(i - j) for j in axes(A, 2), i in axes(A, 1) if !iszero(A[i, j]))
    diagnostics = merge(metadata, (;
        method, requested_order=order, row_orders, positive_norm=true,
        nullity=length(x) - _rank_of_weighted_derivative(Q),
        nnz_H=count_nonzero(H), nnz_Q=count_nonzero(Q), nnz_D=count_nonzero(D),
        bandwidth_H=bandwidth(H), bandwidth_Q=bandwidth(Q), bandwidth_D=bandwidth(D),
        norm_condition=isdiag(H) ? Float64(maximum(diag(H)) / minimum(diag(H))) :
                                  cond(Matrix{Float64}(H)), exact=true,
    ))
    return FunctionSpaceSBPOperator(1, x, D, H, B, Q,
                                    ["x^$k" for k in 0:order], diagnostics)
end

"""
    polynomial_element(order; domain=(0, 1))

Exact Lagrange differentiation and polynomial mass on one uniform element
with `order + 1` nodes. Scaling is by element size.
"""
function polynomial_element(order::Integer; domain::Tuple{<:Real,<:Real}=(0, 1))
    order >= 1 || throw(ArgumentError("order must be positive"))
    a, b = _interval(domain)
    ξ = collect(range(zero(R), one(R); length=order + 1))
    V, V′ = polynomial_matrices(ξ, order, 1)
    C = V \ Matrix{R}(I, order + 1, order + 1)
    G = R[1 // (i + j + 1) for i in 0:order, j in 0:order]
    H = (b - a) .* (transpose(C) * G * C)
    D = (V′ * C) ./ (b - a)
    _positive_definite(H) || error("polynomial element norm failed exact positivity")
    x = a .+ (b - a) .* ξ
    B = _boundary_matrix(length(ξ))
    return _operator(H, D, B, H * D, x, order, :polynomial_element)
end

function _interval(domain)
    a, b = domain
    isfinite(a) && isfinite(b) && b > a ||
        throw(ArgumentError("domain endpoints must be finite and ordered"))
    return _exact(a), _exact(b)
end

function _matrices(system, z)
    k = length(system.norm_orbits)
    Hhat = _function_space_matrix(system.Hfixed, system.norm_orbits, view(z, 1:k), 1)
    Q = _function_space_matrix(system.B / 2, system.skew_orbits, view(z, (k + 1):length(z)), -1)
    return Hhat, Q
end

"""
    stencil_candidate(system, solution=solve_constraints(system.constraints); parameters, domain)

Certify a specific member of the affine experimental system. An indefinite
norm returns `nothing`. The report includes row accuracies and nullity.
"""
function stencil_candidate(system::StencilExperiment,
                           solution::AffineSolution{R}=solve_constraints(system.constraints);
                           parameters=R[], domain::Tuple{<:Real,<:Real}=(0, 1))
    z = evaluate(solution, parameters)
    return _candidate(system, z, _interval(domain))
end

function _candidate(system, z, interval)
    system.constraints.A * z == system.constraints.b ||
        throw(ArgumentError("coefficients do not satisfy the experimental system"))
    Hhat, Q = _matrices(system, z)
    _positive_definite(Hhat) || return nothing
    a, b = interval
    n = size(Hhat, 1)
    h = (b - a) / (n - 1)
    x = [a + h * i for i in 0:(n - 1)]
    H = h .* Hhat
    D = H \ Q
    metadata = (; period=system.period, closure_depth=system.closure_depth,
                stencil_radius=system.stencil_radius, norm=system.norm,
                norm_bandwidth=system.norm_bandwidth,
                closure_norm_bandwidth=system.closure_norm_bandwidth,
                closure_rows=vcat(1:system.closure_depth, (n - system.closure_depth + 1):n),
                spacing=h)
    return _operator(H, D, system.B, Q, x, system.order, :stencil_experiment; metadata)
end

"""
    refine_stencil(system, operator; N, domain=(first(operator.x), last(operator.x)))

Transfer the same experimental coefficients to a larger grid with the same
reflection phases. Recheck exactness, positivity, and the nullspace there.
"""
function refine_stencil(system::StencilExperiment, operator::FunctionSpaceSBPOperator{R};
                        N::Integer, domain::Tuple{<:Real,<:Real}=(first(operator.x), last(operator.x)))
    n = size(system.B, 1)
    length(operator.x) == n || throw(DimensionMismatch("operator must match the original grid"))
    N > n && mod(N - n, system.period) == 0 ||
        throw(ArgumentError("increase N by a positive multiple of the period"))
    refined = assemble_stencil_experiment(system.order; N, period=system.period,
        closure_depth=system.closure_depth, stencil_radius=system.stencil_radius,
        norm=system.norm, norm_bandwidth=system.norm_bandwidth,
        closure_norm_bandwidth=system.closure_norm_bandwidth)
    function transfer(orbits, refined_orbits, matrix)
        m, q = system.closure_depth, system.period
        key(orbit, count) = _orbit_key(first(orbit)..., count, m, q)
        values = Dict(key(orbit, n) => matrix[first(orbit)...] for orbit in orbits)
        length(values) == length(refined_orbits) &&
            all(orbit -> haskey(values, key(orbit, N)), refined_orbits) ||
            throw(ArgumentError("the original grid must sample every transferred phase"))
        return R[values[key(orbit, N)] for orbit in refined_orbits]
    end
    h = (last(operator.x) - first(operator.x)) / (n - 1)
    weights = transfer(system.norm_orbits, refined.norm_orbits, operator.H / h)
    skew = transfer(system.skew_orbits, refined.skew_orbits, operator.Q)
    original = vcat([operator.H[first(orbit)...] / h for orbit in system.norm_orbits],
                    [operator.Q[first(orbit)...] for orbit in system.skew_orbits])
    Hhat, Q = _matrices(system, original)
    Hhat * h == operator.H && Q == operator.Q ||
        throw(ArgumentError("operator does not match the original stencil structure"))
    return _candidate(refined, vcat(weights, skew), _interval(domain))
end

"""
    search_stencil(system; domain=(0, 1), max_directions=32)

Try a bounded set of exact affine candidates and retain a positive
constant-kernel operator with exact order at all rows when one is found.
Returns the best checked positive candidate and the complete affine family
even if the search does not find exact order at every row.
"""
function search_stencil(system::StencilExperiment;
                        domain::Tuple{<:Real,<:Real}=(0, 1), max_directions::Integer=32)
    max_directions >= 0 || throw(ArgumentError("max_directions must be nonnegative"))
    interval = _interval(domain)
    solution = solve_constraints(system.constraints)
    isnothing(solution.particular) &&
        return (; status=:inconsistent_linear_system, operator=nothing, solution, checked=0)

    particular = evaluate(solution)
    target = R[i == j ? 1 : 0 for (i, j) in first.(system.norm_orbits)]
    metric = R[length(orbit) * (i == j ? 1 : 2)
               for (orbit, (i, j)) in zip(system.norm_orbits, first.(system.norm_orbits))]
    near = _near_norm_solution(solution, target, metric, :exact)
    candidates = Iterators.flatten(((particular, near),
        (center + step * v for v in Iterators.take(eachcol(solution.nullspace), max_directions)
         for center in (near, particular) for step in (1//2, -1//2, 1, -1))))
    best, best_score, checked = nothing, (-1, -1), 0
    for z in candidates
        checked += 1
        candidate = _candidate(system, z, interval)
        isnothing(candidate) && continue
        diagnostics = candidate.diagnostics
        score = (Int(diagnostics.nullity == 1), count(==(system.order), diagnostics.row_orders))
        if score == (1, length(candidate.x))
            return (; status=:success, operator=candidate, solution, checked)
        elseif score > best_score
            best, best_score = candidate, score
        end
    end
    status = isnothing(best) ? :positive_norm_not_found :
             best_score[1] == 0 ? :nullspace_not_certified : :exact_order_not_found
    return (; status, operator=best, solution, checked)
end

"""
    wider_stencil_obstruction(order; radius, closure_depth, N)

Leave every centered interior tap unknown and test whether the degree
`order + 1` residual vanishes on the complete affine family.
"""
function wider_stencil_obstruction(order::Integer;
    radius::Integer=order ÷ 2 + 1, closure_depth::Integer=max(radius, order + 1),
    N::Integer=2 * closure_depth + 2 * radius + 1,
)
    system = assemble_stencil_experiment(order; N, closure_depth, stencil_radius=radius,
                                         norm=:identity_interior)
    solution = solve_constraints(system.constraints)
    isnothing(solution.particular) &&
        return (; status=:inconsistent_linear_system, forced_zero=nothing, system, solution)
    X, Y = polynomial_matrices(0:(N - 1), order + 1, 1)
    X, Y = X[:, end:end], Y[:, end:end]
    m = system.closure_depth
    interior = (m + 1):(N - m)
    moments = hcat(-_actions(Y, system.norm_orbits, 1),
                    _actions(X, system.skew_orbits, -1))[interior, :]
    offset = ((system.B / 2) * X - system.Hfixed * Y)[interior]
    forced_zero = all(iszero, moments * solution.particular + offset) &&
                  all(iszero, moments * solution.nullspace)
    return (; status=:consistent, forced_zero, system, solution)
end
