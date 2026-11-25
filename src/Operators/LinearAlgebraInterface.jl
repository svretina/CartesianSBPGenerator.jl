derivative_matrix(operator::SBPOperator) = operator.D
mass_matrix(operator::SBPOperator) = operator.H
boundary_matrix(operator::SBPOperator) = operator.B

function gradient_matrix(operator::SBPOperator)
    operator.derivative_order == 1 || throw(ArgumentError(
        "gradient_matrix is only defined for first-derivative operators",
    ))
    return derivative_matrix(operator)
end

function laplacian_matrix(operator::SBPOperator)
    operator.derivative_order == 2 || throw(ArgumentError(
        "laplacian_matrix is only defined for second-derivative operators",
    ))
    return derivative_matrix(operator)
end
