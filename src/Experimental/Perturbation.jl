function _perturbation_tile(p)
    d = R[(-1)^(p + 1 - i) * binomial(big(p + 1), i) for i in 0:(p + 1)]
    u = vcat(d, zero(R))
    v = vcat(zero(R), d)
    return u*v' - v*u'
end

function _perturbation_bulk_taps(tile, offset)
    return sum((-1)^(i-1) * tile[i,i+offset]
               for i in max(1,1-offset):min(size(tile,1),size(tile,1)-offset))
end

"""Signed dimensionless coefficient of the degree-`order+1` Taylor term."""
function leading_error_coefficients(op, order::Integer)
    order >= 1 || throw(ArgumentError("order must be positive"))
    n = length(op.x)
    h = (last(op.x) - first(op.x)) / (n - 1)
    degree = order + 1
    return [h * sum(op.D[i,j] * big(j-i)^degree for j in 1:n) / factorial(big(degree))
            for i in 1:n]
end

"""
    analytical_perturbation(baseline; N, epsilon, domain=(0,1))
    analytical_perturbation(order; N, epsilon, domain=(0,1))

Add an overlapping alternating skew perturbation to an exact diagonal-norm
SBP(2p,p) baseline. Use odd N to preserve reflection. The default amplitude
is 1/(12*2^(2p)). This is an independent experimental construction.
"""
function analytical_perturbation(
    baseline::SBPFamily{R};
    N::Integer,
    epsilon::Union{Nothing,Integer,Rational}=nothing,
    domain::Tuple{<:Real,<:Real}=(0,1),
)
    p = baseline.achieved_boundary_order
    p >= 2 && iseven(p) && baseline.achieved_interior_order == 2p &&
        baseline.norm === :diagonal ||
        throw(ArgumentError("baseline must be an exact diagonal SBP(2p,p) family with positive even p"))
    isodd(N) || throw(ArgumentError("analytical perturbation requires odd N for reflection"))
    eps = isnothing(epsilon) ? one(R) / (12*big(2)^(2p)) : _exact(epsilon)
    iszero(eps) && throw(ArgumentError("epsilon must be nonzero"))
    depth = max(baseline.closure_depth, p+2)
    N >= 2depth + 2 || throw(ArgumentError("grid must contain a full untruncated bulk period"))
    base = instantiate(baseline; N, domain=_interval(domain))
    H, B, x = base.H, base.B, base.x
    tile = _perturbation_tile(p)
    width = size(tile,1)
    K = zeros(R,N,N)
    for j in 1:(N-width+1)
        rows = j:(j+width-1)
        K[rows,rows] .+= (-1)^(j-1) .* tile
    end
    X, _ = polynomial_matrices(0:N-1,p,1)
    K == -K' && all(iszero,K*X) || error("invalid polynomial-annihilating perturbation")
    Q = weighted_derivative_matrix(base) + eps*K
    D = Diagonal(diag(H)) \ Q
    offsets = -(p+2):(p+2)
    base_taps = coefficients(baseline.interior_stencil)
    bulk_taps = ntuple(2) do phase
        R[(abs(j)<=p ? base_taps[j+p+1] : 0) +
          eps*(-1)^(phase-1)*_perturbation_bulk_taps(tile,j) for j in offsets]
    end
    h = (last(x)-first(x))/(N-1)
    for i in (depth+1):(N-depth)
        h .* D[i, i .+ offsets] == bulk_taps[mod1(i,2)] ||
            error("overlapping tiles do not reproduce the analytical bulk taps")
    end
    D == -reverse(D;dims=(1,2)) || error("analytical perturbation lost reflection")
    all(>(0),diag(H)) || error("analytical baseline lost positivity")
    metadata = (; epsilon=eps, closure_depth=depth, stencil_radius=p+2,
                baseline_closure_depth=baseline.closure_depth, period=2, bulk_taps, offsets,
                spacing=h)
    op = _operator(H,D,B,Q,x,p,:analytical_perturbation; metadata)
    op.diagnostics.nullity == 1 && all(==(p),op.diagnostics.row_orders) ||
        throw(ArgumentError("epsilon does not yield exact requested order and a constant-only nullspace on this grid"))
    return op
end

function analytical_perturbation(order::Integer; N::Integer, kwargs...)
    order >= 2 && iseven(order) || throw(ArgumentError("order must be positive and even"))
    baseline = design_sbp(2order, order; norm=:diagonal, skew_reconstruction=:analytical)
    return analytical_perturbation(baseline; N, kwargs...)
end
