# Keep only independent directions in the projection of the affine family
# onto its norm coordinates. Skew-only directions do not affect feasibility.
function _norm_directions(solution, count, backend)
    directions = solution.nullspace[1:count, :]
    isempty(directions) && return Int[]
    result = solve_constraints(LinearConstraintSystem(
        directions, zeros(eltype(directions), size(directions, 1)),
    ); backend)
    return isnothing(result.particular) ? nothing : result.pivot_columns
end

# Each row represents a ⋅ θ + b > 0. Eliminate the last parameter, then
# pick a rational point in the open interval defined by its remaining bounds.
function _strict_feasibility(A::Matrix{T}, b::Vector{T}, budget) where {T<:Rational}
    m, n = size(A)
    active = Int[]
    for i in 1:m
        if all(iszero, view(A, i, :))
            b[i] > 0 || return :infeasible, T[]
        else
            push!(active, i)
        end
    end
    n == 0 && return :feasible, T[]
    lower = [i for i in active if A[i, n] > 0]
    upper = [i for i in active if A[i, n] < 0]
    independent = [i for i in active if iszero(A[i, n])]
    rows = length(lower) * length(upper) + length(independent)
    budget[] -= rows
    budget[] >= 0 || return :unresolved, T[]
    reduced_A = zeros(T, rows, n - 1)
    reduced_b = zeros(T, rows)
    k = 0
    for i in independent
        k += 1
        reduced_A[k, :] .= view(A, i, 1:(n - 1))
        reduced_b[k] = b[i]
    end
    for i in lower, j in upper
        k += 1
        reduced_A[k, :] .= view(A, i, 1:(n - 1)) ./ A[i, n] .-
                           view(A, j, 1:(n - 1)) ./ A[j, n]
        reduced_b[k] = b[i] / A[i, n] - b[j] / A[j, n]
    end
    status, θ = _strict_feasibility(reduced_A, reduced_b, budget)
    status === :feasible || return status, θ
    bound(i) = -(dot(view(A, i, 1:(n - 1)), θ) + b[i]) / A[i, n]
    lo = isempty(lower) ? nothing : maximum(bound, lower)
    hi = isempty(upper) ? nothing : minimum(bound, upper)
    value = if isnothing(lo)
        isnothing(hi) ? zero(T) : hi - one(T)
    elseif isnothing(hi)
        lo + one(T)
    else
        (lo + hi) / 2
    end
    return :feasible, vcat(θ, value)
end

function _select_positive_norm(system::FirstDerivativeSystem{T}, solution, margin) where {T}
    threshold = T <: Rational ? _exact(margin) : T(margin)
    function acceptable(z)
        block = _left_norm(system, z)
        return _positive_definite(block - threshold * I)
    end
    particular = evaluate(solution)
    acceptable(particular) && return :success, particular
    entries = system.norm_entries
    directions = _norm_directions(solution, length(entries), _backend(system))
    if T <: Rational
        isempty(directions) && return :certified_positivity_infeasible, nothing
        for (k, (i, j)) in enumerate(entries)
            if i == j && particular[k] <= threshold &&
               all(iszero, view(solution.nullspace, k, :))
                return :certified_positivity_infeasible, nothing
            end
        end
    end
    reference = _near_identity_solution(system, solution)
    acceptable(reference) && return :success, reference
    isnothing(directions) && return :numerically_uncertain, nothing

    if T <: Rational && system.norm === :diagonal
        status, θ = _strict_feasibility(
            solution.nullspace[1:length(entries), directions],
            particular[1:length(entries)] .- threshold, Ref(50_000),
        )
        status === :infeasible && return :certified_positivity_infeasible, nothing
        if status === :feasible
            candidate = particular + solution.nullspace[:, directions] * θ
            acceptable(candidate) || error("exact positivity feasibility certificate failed")
            return :success, candidate
        end
    else
        for j in Iterators.take(directions, 64), step in (1//2, -1//2, 1, -1, 2, -2)
            candidate = reference + T(step) .* view(solution.nullspace, :, j)
            acceptable(candidate) && return :success, candidate
        end
    end
    return :positive_norm_not_found, nothing
end
