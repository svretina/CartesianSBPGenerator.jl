function instantiate(
    family::SBPFamily{T};
    N::Integer,
    domain::Tuple{<:Real,<:Real}=(0, 1),
) where {T}
    N >= 2 * family.closure_width + 1 || throw(ArgumentError(
        "N must be at least $(2 * family.closure_width + 1) for this family",
    ))
    left, right = domain
    right > left || throw(ArgumentError("domain must satisfy right > left"))
    S = left isa AbstractFloat || right isa AbstractFloat ?
        promote_type(typeof(left), typeof(right)) : T
    converted_left, converted_right = convert(S, left), convert(S, right)
    h = (converted_right - converted_left) / (N - 1)
    x = [converted_left + (i - 1) * h for i in 1:N]
    m, w = family.closure_depth, family.closure_width
    left_rows, right_rows = 1:m, (N - m + 1):N
    left_block = family.norm === :diagonal ?
        Diagonal(S.(family.left_weights)) : S.(family.left_norm)
    right_block = reverse(reverse(Matrix(left_block); dims=1); dims=2)
    H = Matrix{S}(I, N, N) .* h
    H[left_rows, left_rows] .= h .* left_block
    H[right_rows, right_rows] .= h .* right_block

    Q = zeros(S, N, N)
    c = S.(coefficients(family.interior_stencil))
    radius = family.interior_stencil.radius
    for row in (m + 1):(N - m)
        Q[row, (row - radius):(row + radius)] .= c
    end
    Q[1:m, 1:w] .= S.(family.left_Q)
    Q[(N - m + 1):N, (N - w + 1):N] .=
        -S.(reverse(reverse(family.left_Q; dims=1); dims=2))
    # Only the two boundary blocks require a solve.
    D = Q / h
    D[left_rows, 1:w] .= left_block \ D[left_rows, 1:w]
    D[right_rows, (N - w + 1):N] .= right_block \ D[right_rows, (N - w + 1):N]
    B = zeros(S, N, N)
    B[1, 1] = -one(S)
    B[N, N] = one(S)
    diagnostics = merge(family.diagnostics, (N=N, spacing=h,))
    return SBPOperator(
        family.derivative_order,
        family.requested_interior_order,
        family.requested_boundary_order,
        family.achieved_interior_order,
        family.achieved_boundary_order,
        x, Matrix(D), H, B, (Q=Q,), diagnostics,
    )
end

function generate_sbp(
    interior_order::Integer,
    boundary_order::Integer;
    N::Integer=101,
    domain::Tuple{<:Real,<:Real}=(0, 1),
    kwargs...,
)
    family = design_sbp(interior_order, boundary_order; kwargs...)
    return instantiate(family; N=N, domain=domain)
end
