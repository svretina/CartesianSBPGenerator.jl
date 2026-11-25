struct FirstDerivativeSystem{T}
    constraints::LinearConstraintSystem{T,Matrix{T},Vector{T}}
    spec::SBPSpec
    closure_depth::Int
    closure_width::Int
    canonical_size::Int
    interior_stencil::CenteredStencil{T}
    norm::Symbol
    norm_entries::Vector{Tuple{Int,Int}}
end

function _default_closure_depth(boundary_order::Int)
    return boundary_order == 1 ? 1 : 2 * boundary_order
end

function _first_derivative_parameterizations(
    stencil::CenteredStencil{T}, closure_depth::Int, closure_width::Int, n::Int,
    norm_entries,
) where {T}
    number_of_h_variables = length(norm_entries)
    number_of_q_variables = closure_depth * closure_width
    number_of_variables = number_of_h_variables + number_of_q_variables
    H0 = Matrix{T}(I, n, n)
    for i in 1:closure_depth
        H0[i, i] = zero(T)
        H0[n + 1 - i, n + 1 - i] = zero(T)
    end
    Hbasis = [zeros(T, n, n) for _ in 1:number_of_variables]
    for (variable, (i, j)) in enumerate(norm_entries)
        Hbasis[variable][i, j] = Hbasis[variable][j, i] = one(T)
        Hbasis[variable][n + 1 - i, n + 1 - j] =
            Hbasis[variable][n + 1 - j, n + 1 - i] = one(T)
    end

    Q0 = zeros(T, n, n)
    c = coefficients(stencil)
    radius = stencil.radius
    for row in (closure_depth + 1):(n - closure_depth)
        Q0[row, (row - radius):(row + radius)] .= c
    end
    Qbasis = [zeros(T, n, n) for _ in 1:number_of_variables]
    variable = number_of_h_variables
    for column in 1:closure_width, row in 1:closure_depth
        variable += 1
        Qbasis[variable][row, column] = one(T)
        Qbasis[variable][n + 1 - row, n + 1 - column] = -one(T)
    end
    names = vcat(
        [i == j ? Symbol(:w_, i - 1) : Symbol(:h_, i - 1, :_, j - 1)
         for (i, j) in norm_entries],
        [Symbol(:q_, row - 1, :_, column - 1)
         for column in 1:closure_width for row in 1:closure_depth],
    )
    return H0, Hbasis, Q0, Qbasis, names
end

function _append_equation!(rows, right_hand_sides, coefficients, right_hand_side)
    push!(rows, coefficients)
    push!(right_hand_sides, right_hand_side)
    return nothing
end

"""Assemble the exact linear equations for a first-derivative SBP closure."""
function assemble_constraints(spec::SBPSpec)
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
    T = eltype(stencil.particular)
    norm_entries = [(i, i) for i in 1:closure_depth]
    if norm === :restricted_full
        bandwidth = something(spec.norm_bandwidth, closure_depth - 1)
        append!(norm_entries, [
            (i, j) for j in 3:closure_depth for i in max(2, j - bandwidth):(j - 1)
        ])
    end
    H0, Hbasis, Q0, Qbasis, names = _first_derivative_parameterizations(
        stencil, closure_depth, closure_width, n, norm_entries,
    )
    number_of_variables = length(names)
    rows = Vector{T}[]
    right_hand_sides = T[]
    B = zeros(T, n, n)
    B[1, 1] = -one(T)
    B[n, n] = one(T)

    # Q + Q' = B. Only retain nontrivial equations.
    for column in 1:n, row in 1:column
        constant_term = Q0[row, column] + Q0[column, row]
        equation = T[
            Qbasis[v][row, column] + Qbasis[v][column, row]
            for v in 1:number_of_variables
        ]
        rhs = B[row, column] - constant_term
        (any(!iszero, equation) || !iszero(rhs)) &&
            _append_equation!(rows, right_hand_sides, equation, rhs)
    end

    # Boundary exactness Q*X = H*X' through degree p_b.
    nodes = collect(0:(n - 1))
    X, Xd = polynomial_matrices(nodes, spec.boundary_order, 1; arithmetic=T)
    for degree in 0:spec.boundary_order, row in 1:closure_depth
        constant_term = dot(Q0[row, :], X[:, degree + 1]) -
                        dot(H0[row, :], Xd[:, degree + 1])
        equation = T[
            dot(Qbasis[v][row, :], X[:, degree + 1]) -
            dot(Hbasis[v][row, :], Xd[:, degree + 1])
            for v in 1:number_of_variables
        ]
        _append_equation!(rows, right_hand_sides, equation, -constant_term)
    end

    A = isempty(rows) ? zeros(T, 0, number_of_variables) : permutedims(hcat(rows...))
    constraints = LinearConstraintSystem(A, right_hand_sides, names)
    return FirstDerivativeSystem(
        constraints, spec, closure_depth, closure_width, n, stencil, norm, norm_entries,
    )
end

solve_constraints(system::FirstDerivativeSystem; kwargs...) =
    solve_constraints(system.constraints; kwargs...)
