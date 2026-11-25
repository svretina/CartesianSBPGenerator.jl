function _falling_factorial(degree::Integer, derivative_order::Integer)
    derivative_order > degree && return BigInt(0)
    result = BigInt(1)
    for value in (degree - derivative_order + 1):degree
        result *= value
    end
    return result
end

"""Evaluate the `d`th derivative of `x^degree` at `x`."""
function derivative_monomial(
    x,
    degree::Integer,
    derivative_order::Integer;
    arithmetic::Type=promote_type(typeof(x), BigRational),
)
    degree >= 0 || throw(ArgumentError("degree must be nonnegative"))
    derivative_order >= 0 || throw(ArgumentError("derivative_order must be nonnegative"))
    derivative_order > degree && return zero(arithmetic)
    factor = convert(arithmetic, _falling_factorial(degree, derivative_order))
    return factor * convert(arithmetic, x)^(degree - derivative_order)
end

"""Construct polynomial value and derivative matrices on a grid."""
function polynomial_matrices(
    nodes::AbstractVector,
    maximum_degree::Integer,
    derivative_order::Integer;
    arithmetic::Type=BigRational,
)
    maximum_degree >= 0 || throw(ArgumentError("maximum_degree must be nonnegative"))
    derivative_order >= 0 || throw(ArgumentError("derivative_order must be nonnegative"))
    X = zeros(arithmetic, length(nodes), maximum_degree + 1)
    Xd = similar(X)
    for degree in 0:maximum_degree, (row, node) in enumerate(nodes)
        converted_node = convert(arithmetic, node)
        X[row, degree + 1] = converted_node^degree
        Xd[row, degree + 1] = derivative_monomial(
            converted_node, degree, derivative_order; arithmetic=arithmetic,
        )
    end
    return X, Xd
end
