using CartesianSBPGenerator
using CartesianSBPGenerator: Experimental

for order in (2, 4, 6), extra_radius in (1, 2)
    radius = order ÷ 2 + extra_radius
    depth = order + extra_radius
    result = Experimental.wider_stencil_obstruction(order; radius, closure_depth=depth)
    println((; experiment=:wider_centered, order, radius, depth,
             result.status, result.forced_zero, free=nullity(result.solution)))
end

for order in (2, 4, 6)
    operator = Experimental.polynomial_element(order)
    println((; experiment=:polynomial_element, operator.diagnostics...))
end

for (order, period, norm_bandwidth, depth, radius, N) in (
    (2, 2, 0, 5, 2, 17), (2, 3, 0, 5, 3, 19), (2, 4, 0, 5, 2, 23),
    (4, 2, 0, 7, 4, 21), (4, 3, 0, 7, 4, 24), (4, 4, 0, 7, 4, 27),
    (4, 2, 0, 9, 4, 25), (4, 4, 0, 9, 4, 31),
    (6, 2, 0, 7, 4, 20), (6, 3, 0, 7, 4, 23), (6, 4, 0, 7, 4, 26),
    (6, 2, 0, 15, 6, 37), (6, 4, 0, 15, 6, 43),
    (2, 2, 1, 3, 2, 12), (4, 2, 1, 5, 3, 16), (6, 2, 1, 7, 4, 20),
    (2, 1, 1, 3, 2, 9), (4, 1, 1, 5, 3, 13), (6, 1, 1, 7, 4, 17),
    (4, 1, 2, 6, 3, 17), (6, 1, 2, 10, 4, 25), (6, 1, 3, 10, 4, 25),
)
    system = Experimental.assemble_stencil_experiment(order; N, period, closure_depth=depth,
                stencil_radius=radius, norm_bandwidth)
    seconds = @elapsed result = Experimental.search_stencil(system; max_directions=8)
    println((; experiment=:stencil, order, period, norm_bandwidth, depth, radius, N,
             result.status, result.checked, seconds))
end

for (order, period, norm_bandwidth, depth, radius, N) in (
    (2, 2, 0, 5, 2, 17), (4, 2, 0, 9, 4, 25), (6, 2, 0, 15, 6, 37),
    (4, 1, 2, 6, 3, 17), (6, 1, 3, 10, 4, 25),
)
    system = Experimental.assemble_stencil_experiment(order; N, period, closure_depth=depth,
                stencil_radius=radius, norm_bandwidth)
    result = Experimental.search_stencil(system)
    result.status === :success || continue
    for n in (N, 2N - 1, 4N - 3)
        seconds = @elapsed op = n == N ? result.operator : Experimental.refine_stencil(system, result.operator; N=n)
        isnothing(op) && error("refined norm failed positivity")
        x = Float64.(op.x)
        data = exp.(x)
        residual = abs.(Float64.(op.D) * data - data)
        d = op.diagnostics
        println((; experiment=:refinement, order, period, N=n, seconds,
                 row_order=extrema(d.row_orders), d.nullity, d.norm_condition,
                 d.nnz_H, d.nnz_Q, d.nnz_D,
                 boundary_error=maximum(residual[d.closure_rows]),
                 interior_error=maximum(residual[(depth + 1):(n - depth)])))
    end
end
