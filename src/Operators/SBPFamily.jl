function _left_norm(system::FirstDerivativeSystem, variables::Vector{T}) where {T}
    block = zeros(T, system.closure_depth, system.closure_depth)
    for (value, (i, j)) in zip(variables, system.norm_entries)
        block[i, j] = block[j, i] = value
    end
    return block
end

# Unpivoted LDLᵀ certifies positivity using rational Schur complements.
function _positive_definite(block::Matrix{<:Rational})
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

_positive_definite(block::Matrix{<:AbstractFloat}) = isposdef(Symmetric(block))

function _near_norm_solution(solution, target, weights, backend)
    # The exact normal equations minimize the Frobenius distance to identity.
    # A singular Gram matrix just means some free directions leave H unchanged.
    directions = solution.nullspace[1:length(target), :]
    residual = target - solution.particular[1:length(target)]
    weighted_transpose = transpose(weights .* directions)
    projection = solve_constraints(LinearConstraintSystem(
        weighted_transpose * directions, weighted_transpose * residual,
    ); backend)
    isnothing(projection.particular) && return evaluate(solution)
    return evaluate(solution, evaluate(projection))
end

function _near_identity_solution(system::FirstDerivativeSystem{T}, solution) where {T}
    target = T[i == j ? 1 : 0 for (i, j) in system.norm_entries]
    weights = T[i == j ? 1 : 2 for (i, j) in system.norm_entries]
    return _near_norm_solution(solution, target, weights, _backend(system))
end

"""
    design_sbp(interior_order, boundary_order; closure_depth=:auto, max_closure_depth=32, ...)

Search increasing closure depths at the requested orders. An explicit depth
tests one size. The search report is available as `family.diagnostics.attempts`
or `SBPConstructionError.attempts`.
"""
function design_sbp(
    interior_order::Integer,
    boundary_order::Integer;
    derivative_order::Integer=1,
    norm::Symbol=:auto,
    order_policy::Symbol=:exact,
    interior_radius::Union{Nothing,Integer,Symbol}=nothing,
    closure_depth::Union{Nothing,Integer,Symbol}=:auto,
    max_closure_depth::Integer=32,
    norm_bandwidth::Union{Nothing,Integer}=nothing,
    positivity_margin::Real=0,
    boundary_derivative_order::Union{Nothing,Integer}=nothing,
    compatibility::Symbol=:independent,
    extra_accuracy::Symbol=:auto,
    symmetric::Bool=true,
    construction::Union{Nothing,Symbol}=nothing,
    method::Union{Nothing,Symbol}=nothing,
    arithmetic::Symbol=:exact,
    basis::Union{Nothing,BasisSamples,PolynomialBasis,MonomialBasis}=nothing,
    skew_reconstruction::Symbol=:linear,
)
    interior_radius === :minimal && (interior_radius = nothing)
    interior_radius isa Symbol && throw(ArgumentError("interior_radius must be an integer or :minimal"))
    closure_depth === :auto && (closure_depth = nothing)
    closure_depth isa Symbol && throw(ArgumentError("closure_depth must be an integer or :auto"))
    options = (; interior_order, boundary_order, derivative_order, order_policy,
               interior_radius, norm_bandwidth, boundary_derivative_order, compatibility,
               extra_accuracy, symmetric)
    spec = _validate_spec(SBPSpec(; options..., norm, closure_depth))
    _strategy(spec) isa FirstDerivativeSBP ||
        throw(ArgumentError("second-derivative constraint assembly is not implemented yet"))
    selected = _closure_method(method, construction)
    selected in (:direct, :reduced) || throw(ArgumentError("construction must be :reduced or :direct"))
    skew_reconstruction in (:linear, :analytical) ||
        throw(ArgumentError("skew_reconstruction must be :linear or :analytical"))
    arithmetic in (:exact, :numerical) || throw(ArgumentError("arithmetic must be :exact or :numerical"))
    isfinite(positivity_margin) && positivity_margin >= 0 ||
        throw(ArgumentError("positivity_margin must be finite and nonnegative"))
    max_closure_depth > 0 || throw(ArgumentError("max_closure_depth must be positive"))

    candidate_order = interior_order
    if order_policy === :allow_upgrade
        minimum_order = norm === :diagonal ? 2 * boundary_order : 2 * (boundary_order ÷ 2 + 1)
        candidate_order = max(interior_order, minimum_order)
    elseif norm === :diagonal && boundary_order > interior_order ÷ 2
        throw(ArgumentError("a first-derivative diagonal norm requires interior_order >= 2*boundary_order"))
    elseif boundary_order >= interior_order
        throw(SBPConstructionError(interior_order, boundary_order, :theoretical_obstruction,
                                   SBPClosureAttempt[]))
    end

    minimum_radius = minimum_centered_radius(1, candidate_order)
    radius = something(interior_radius, minimum_radius)
    radius == minimum_radius ||
        throw(ArgumentError("first-derivative SBP construction requires the minimum interior radius $minimum_radius"))
    if basis isa BasisSamples && isnothing(closure_depth)
        closure_depth = size(basis.X, 1) - radius
    end
    isnothing(closure_depth) || closure_depth >= radius ||
        throw(ArgumentError("closure_depth must be at least the interior radius"))
    if isnothing(closure_depth) && max_closure_depth < radius
        throw(ArgumentError("max_closure_depth must be at least the interior radius"))
    end

    diagonal_orders = boundary_order <= candidate_order ÷ 2
    norms = if norm === :auto
        diagonal_orders && isnothing(norm_bandwidth) ? (:diagonal, :restricted_full) : (:restricted_full,)
    else
        (norm,)
    end
    depths = isnothing(closure_depth) ? (radius:max_closure_depth) : (closure_depth:closure_depth)
    attempts = SBPClosureAttempt[]
    candidate_options = merge(options, (interior_order=candidate_order,))
    for candidate_norm in norms, m in depths
        candidate = SBPSpec(; candidate_options..., norm=candidate_norm, closure_depth=m)
        system = assemble_constraints(candidate; construction=selected, arithmetic, basis)
        family = _design_candidate(system, skew_reconstruction, positivity_margin, attempts,
                                   interior_order)
        isnothing(family) || return family
    end
    status = isnothing(closure_depth) ? :search_limit_reached : last(attempts).status
    throw(SBPConstructionError(interior_order, boundary_order, status, attempts))
end

function _design_candidate(system, skew_reconstruction, positivity_margin, attempts, requested_order)
    solution = solve_constraints(system)
    m, w = system.closure_depth, system.closure_width
    construction = system.construction
    radius = system.interior_stencil.radius
    record(status) = push!(attempts, SBPClosureAttempt(system.norm, m, radius, status))
    if isnothing(solution.particular)
        record(solution.status === :inconsistent ? :inconsistent_linear_system : solution.status)
        return nothing
    end
    status, z = _select_positive_norm(system, solution, positivity_margin)
    if status !== :success
        record(status)
        return nothing
    end
    left_norm = _left_norm(system, z)
    rank, free = solution.rank, nullity(solution)
    norm_rank, norm_nullity, skew_nullity = nothing, nothing, nothing
    if construction === :reduced
        skew = reconstruct_skew(system, left_norm; method=skew_reconstruction)
        if isnothing(skew.particular)
            record(skew.status === :inconsistent ? :validation_failed : skew.status)
            return nothing
        end
        skew_entries = evaluate(skew)
        norm_rank, norm_nullity, skew_nullity = rank, free, nullity(skew)
        rank += skew.rank
        free += skew_nullity
    else
        skew_entries = z[(length(system.norm_entries) + 1):end]
    end
    A = _skew_block(skew_entries, m, eltype(z))
    if !_closure_residual_ok(A * system.X, left_norm * system.Y + system.K)
        record(:validation_failed)
        return nothing
    end
    A[1, 1] = -one(eltype(A)) / 2
    left_Q = hcat(A, system.C)
    record(:success)
    diagnostics = (
        constraint_rank=rank,
        constraint_nullity=free,
        norm_rank=norm_rank,
        norm_nullity=norm_nullity,
        skew_nullity=skew_nullity,
        canonical_size=system.canonical_size,
        construction=construction,
        method=construction === :direct ? :classical : :reduced,
        exact=_backend(system) === :exact,
        attempts=copy(attempts),
    )
    p, b = system.spec.interior_order, system.spec.boundary_order
    return SBPFamily(1, requested_order, b, p, b, system.norm, m, w, diag(left_norm),
                     left_norm, left_Q, system.interior_stencil, diagnostics)
end
