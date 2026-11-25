function _numerical_solve(system::LinearConstraintSystem, atol::Real)
    isfinite(atol) && atol >= 0 || throw(ArgumentError("atol must be finite and nonnegative"))
    A = Matrix{Float64}(system.A)
    b = Vector{Float64}(system.b)
    all(isfinite, A) && all(isfinite, b) ||
        throw(ArgumentError("numerical constraints must be finite in Float64"))
    m, n = size(A)
    particular = zeros(n)
    pivots = Int[]
    if min(m, n) == 0
        basis = Matrix{Float64}(I, n, n)
        rank_A = 0
        matrix_scale = 0.0
    else
        factors = svd(A; full=true)
        matrix_scale = first(factors.S)
        rank_tolerance = max(Float64(atol), max(m, n) * eps(Float64) * matrix_scale)
        rank_A = count(>(rank_tolerance), factors.S)
        if any(s -> rank_tolerance / 2 < s <= 2 * rank_tolerance, factors.S)
            return AffineSolution(:numerically_uncertain, rank_A, pivots, nothing, zeros(n, 0))
        end
        basis = factors.V[:, (rank_A + 1):n]
        if rank_A > 0
            pivoted_qr = qr(A, ColumnNorm())
            pivots = pivoted_qr.p[1:rank_A]
            problem = LinearSolve.LinearProblem(A[:, pivots], b)
            result = LinearSolve.solve(problem, LinearSolve.QRFactorization(ColumnNorm()))
            particular[pivots] = result.u
        end
    end

    scale = matrix_scale * norm(particular) + norm(b)
    tolerance = max(Float64(atol), max(m, n) * eps(Float64) * scale)
    residual = norm(A * particular - b)
    if !isfinite(residual) || !isfinite(tolerance)
        return AffineSolution(:numerically_uncertain, rank_A, pivots, nothing, zeros(n, 0))
    elseif residual > tolerance
        status = residual > 2 * tolerance ? :inconsistent : :numerically_uncertain
        return AffineSolution(status, rank_A, pivots, nothing, zeros(n, 0))
    end
    status = rank_A == n ? :unique : :underdetermined
    return AffineSolution(status, rank_A, pivots, particular, basis)
end
