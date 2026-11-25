"""Minimum radius of a centered stencil for derivative `d` and even accuracy `p`."""
function minimum_centered_radius(derivative_order::Integer, accuracy_order::Integer)
    derivative_order > 0 || throw(ArgumentError("derivative_order must be positive"))
    accuracy_order > 0 || throw(ArgumentError("accuracy_order must be positive"))
    iseven(accuracy_order) || throw(ArgumentError(
        "standard centered stencils require an even accuracy_order",
    ))
    return fld(derivative_order + 1, 2) + fld(accuracy_order, 2) - 1
end

"""Assemble moment and parity equations for a centered finite-difference stencil."""
function moment_system(
    derivative_order::Integer,
    accuracy_order::Integer,
    radius::Integer;
    arithmetic::Type=BigRational,
)
    minimum = minimum_centered_radius(derivative_order, accuracy_order)
    radius >= minimum || throw(ArgumentError(
        "radius $radius is too small; minimum is $minimum",
    ))
    offsets = -radius:radius
    maximum_degree = derivative_order + accuracy_order - 1
    parity_rows = radius + (isodd(derivative_order) ? 1 : 0)
    A = zeros(arithmetic, maximum_degree + 1 + parity_rows, length(offsets))
    b = zeros(arithmetic, size(A, 1))

    for degree in 0:maximum_degree, (column, offset) in enumerate(offsets)
        A[degree + 1, column] = convert(arithmetic, BigInt(offset)^degree)
    end
    b[derivative_order + 1] = convert(arithmetic, factorial(big(derivative_order)))

    row = maximum_degree + 2
    sign = iseven(derivative_order) ? one(arithmetic) : -one(arithmetic)
    for j in 1:radius
        A[row, radius + 1 - j] = one(arithmetic)
        A[row, radius + 1 + j] = -sign
        row += 1
    end
    if isodd(derivative_order)
        A[row, radius + 1] = one(arithmetic)
    end
    names = [Symbol(:c_, offset < 0 ? :m : :p, abs(offset)) for offset in offsets]
    return LinearConstraintSystem(A, b, names)
end
