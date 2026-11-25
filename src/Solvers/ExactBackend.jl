const BigRational = Rational{BigInt}

_exact(x::Rational) = BigRational(BigInt(numerator(x)), BigInt(denominator(x)))
_exact(x::Integer) = BigRational(BigInt(x), BigInt(1))
_exact(x::AbstractFloat) = rationalize(BigInt, x; tol=0)

function _rref_solve(A::Matrix{T}, b::Vector{T}) where {T}
    m, n = size(A)
    augmented = hcat(copy(A), copy(b))
    pivots = Int[]
    pivot_row = 1

    for column in 1:n
        candidate = findfirst(
            row -> !iszero(augmented[row, column]),
            pivot_row:m,
        )
        isnothing(candidate) && continue
        row = pivot_row + candidate - 1
        if row != pivot_row
            augmented[pivot_row, :], augmented[row, :] =
                copy(augmented[row, :]), copy(augmented[pivot_row, :])
        end

        pivot = augmented[pivot_row, column]
        augmented[pivot_row, :] ./= pivot
        for other_row in 1:m
            other_row == pivot_row && continue
            multiplier = augmented[other_row, column]
            iszero(multiplier) && continue
            augmented[other_row, :] .-= multiplier .* augmented[pivot_row, :]
        end
        push!(pivots, column)
        pivot_row += 1
        pivot_row > m && break
    end

    for row in 1:m
        zero_left = all(column -> iszero(augmented[row, column]), 1:n)
        if zero_left && !iszero(augmented[row, n + 1])
            return AffineSolution{T}(:inconsistent, length(pivots), pivots,
                                     nothing, zeros(T, n, 0))
        end
    end

    particular = zeros(T, n)
    for (row, column) in enumerate(pivots)
        particular[column] = augmented[row, n + 1]
    end
    free_columns = setdiff(collect(1:n), pivots)
    basis = zeros(T, n, length(free_columns))
    for (basis_column, free_column) in enumerate(free_columns)
        basis[free_column, basis_column] = one(T)
        for (row, pivot_column) in enumerate(pivots)
            basis[pivot_column, basis_column] = -augmented[row, free_column]
        end
    end
    status = isempty(free_columns) ? :unique : :underdetermined
    return AffineSolution{T}(status, length(pivots), pivots, particular, basis)
end

"""Solve `A*z=b`, retaining the full affine solution family.

The `:exact` backend promotes input to `Rational{BigInt}`. The `:numerical`
backend uses LinearSolve in Float64 and checks the rank and residual independently.
"""
function solve_constraints(
    system::LinearConstraintSystem;
    backend::Symbol=:exact,
    atol::Real=1.0e-12,
)
    backend in (:exact, :numerical) ||
        throw(ArgumentError("backend must be :exact or :numerical"))
    backend === :numerical && return _numerical_solve(system, atol)
    A = Matrix{BigRational}(_exact.(system.A))
    b = Vector{BigRational}(_exact.(system.b))
    return _rref_solve(A, b)
end
