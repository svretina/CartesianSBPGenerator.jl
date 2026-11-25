function _left_norm(system::FirstDerivativeSystem, variables::Vector{T}) where {T}
    block = zeros(T, system.closure_depth, system.closure_depth)
    for (value, (i, j)) in zip(variables, system.norm_entries)
        block[i, j] = block[j, i] = value
    end
    return block
end

# Unpivoted LDLᵀ certifies positivity using rational Schur complements.
function _exact_positive_definite(block::Matrix{<:Rational})
    remainder = copy(block)
    for i in axes(block, 1)
        pivot = remainder[i, i]
        pivot > 0 || return false
        for j in (i + 1):size(block, 1), k in (i + 1):size(block, 2)
            remainder[j, k] -= remainder[j, i] * remainder[i, k] / pivot
        end
    end
    return true
end

function _near_identity_solution(system::FirstDerivativeSystem{T}, solution) where {T}
    # The exact normal equations minimize the Frobenius distance to identity.
    # A singular Gram matrix just means some free directions leave H unchanged.
    entries = system.norm_entries
    directions = solution.nullspace[1:length(entries), :]
    weights = T[i == j ? 1 : 2 for (i, j) in entries]
    target = T[i == j ? 1 : 0 for (i, j) in entries]
    residual = target - solution.particular[1:length(entries)]
    weighted_transpose = transpose(weights .* directions)
    projection = solve_constraints(LinearConstraintSystem(
        weighted_transpose * directions, weighted_transpose * residual,
    ))
    return evaluate(solution, evaluate(projection))
end

function design_sbp(
    interior_order::Integer,
    boundary_order::Integer;
    derivative_order::Integer=1,
    norm::Symbol=:auto,
    order_policy::Symbol=:exact,
    interior_radius::Union{Nothing,Integer}=nothing,
    closure_depth::Union{Nothing,Integer}=nothing,
    norm_bandwidth::Union{Nothing,Integer}=nothing,
    boundary_derivative_order::Union{Nothing,Integer}=nothing,
    compatibility::Symbol=:independent,
    extra_accuracy::Symbol=:auto,
    symmetric::Bool=true,
)
    spec = SBPSpec(
        interior_order=interior_order,
        boundary_order=boundary_order,
        derivative_order=derivative_order,
        norm=norm,
        order_policy=order_policy,
        interior_radius=interior_radius,
        closure_depth=closure_depth,
        norm_bandwidth=norm_bandwidth,
        boundary_derivative_order=boundary_derivative_order,
        compatibility=compatibility,
        extra_accuracy=extra_accuracy,
        symmetric=symmetric,
    )
    system = assemble_constraints(spec)
    solution = solve_constraints(system; backend=:exact)
    solution.status === :inconsistent && error(
        "the $(system.norm) first-derivative equations are inconsistent " *
        "at closure depth $(system.closure_depth)",
    )
    z = evaluate(solution)
    m, w = system.closure_depth, system.closure_width
    left_norm = _left_norm(system, z)
    if !_exact_positive_definite(left_norm) && nullity(solution) > 0
        z = _near_identity_solution(system, solution)
        left_norm = _left_norm(system, z)
    end
    _exact_positive_definite(left_norm) || error(
        "norm selection found no positive-definite candidate for $(system.norm) " *
        "at closure depth $m",
    )
    left_Q = reshape(z[(length(system.norm_entries) + 1):end], m, w)
    diagnostics = (
        constraint_rank=solution.rank,
        constraint_nullity=nullity(solution),
        canonical_size=system.canonical_size,
        exact=true,
    )
    return SBPFamily(
        1, interior_order, boundary_order, interior_order, boundary_order,
        system.norm, m, w, diag(left_norm), left_norm, left_Q, system.interior_stencil, diagnostics,
    )
end
