"""
    alternating_characteristics(system, operator)

Extract exact dimensionless diagonal weights and stencils in grid-index
phases A (even zero-based index) and B (odd zero-based index). Check that
the extracted stencils reproduce every interior row.
"""
function alternating_characteristics(system::StencilExperiment, op::FunctionSpaceSBPOperator{R})
    system.period == 2 && isdiag(op.H) ||
        throw(ArgumentError("analysis requires a period-two diagonal-norm stencil"))
    n, m, s = size(system.B, 1), system.closure_depth, system.stencil_radius
    length(op.x) == n || throw(DimensionMismatch("operator must match the system grid"))
    h = (last(op.x) - first(op.x)) / (n - 1)
    bulk = (m + 1):(n - m)
    phase(i) = mod1(i, 2)
    representatives = ntuple(r -> first(i for i in bulk if phase(i) == r), 2)
    offsets = -s:s
    weights = ntuple(r -> op.H[representatives[r], representatives[r]] / h, 2)
    stencils = ntuple(r -> [h * op.D[representatives[r], representatives[r] + j] for j in offsets], 2)
    for i in bulk
        r = phase(i)
        op.H[i, i] / h == weights[r] && h .* op.D[i, (i - s):(i + s)] == stencils[r] &&
            count(!iszero, op.D[i, :]) == count(!iszero, stencils[r]) ||
            error("operator contains a nonperiodic interior row")
    end
    p = system.order
    moments = ntuple(r -> [sum(stencils[r] .* big.(offsets) .^ k) for k in 0:(p + 1)], 2)
    leading = ntuple(r -> moments[r][end] / factorial(big(p + 1)), 2)
    wmoment = sum(weights[r] * moments[r][end] for r in 1:2)
    weighted_skew = all(weights[r] * stencils[r][j + s + 1] ==
                       -weights[mod1(r + j, 2)] * stencils[mod1(r + j, 2)][s + 1 - j]
                       for r in 1:2, j in offsets)
    cut_flux = [sum(op.Q[i, j] * big(j - i)^p for i in 1:k, j in (k + 1):n) for k in bulk]
    parity_compatible = all(c == -reverse(c) for c in stencils)
    closure = vcat(1:m, (n - m + 1):n)
    boundary_width = maximum(maximum(findall(!iszero, op.D[i, :])) -
                             minimum(findall(!iszero, op.D[i, :])) + 1 for i in closure)
    row_nnz = [count(!iszero, op.D[i, :]) for i in 1:n]
    return (; order=p, offsets, weights, stencils, moments, leading, weighted_moment=wmoment,
            weighted_skew, cut_flux, parity_compatible, boundary_width, row_nnz,
            maximum_coefficient=maximum(abs, op.D) * h)
end

"""
    bloch_symbol(data, theta; derivative=false)

Dimensionless two-phase symbol in the node gauge, with offsets carrying
the full phase factor. The derivative option differentiates in `theta`.
"""
function bloch_symbol(data, theta::T; derivative::Bool=false) where {T<:AbstractFloat}
    matrix = zeros(Complex{T}, 2, 2)
    for r in 1:2, (j, weight) in zip(data.offsets, data.stencils[r])
        term = T(weight) * cis(j * theta)
        matrix[r, mod1(r + j, 2)] += derivative ? im * j * term : term
    end
    return matrix
end
