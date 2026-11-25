derivative_matrix(operator::Union{SBPOperator,FunctionSpaceSBPOperator}) = operator.D
mass_matrix(operator::Union{SBPOperator,FunctionSpaceSBPOperator}) = operator.H
boundary_matrix(operator::Union{SBPOperator,FunctionSpaceSBPOperator}) = operator.B
weighted_derivative_matrix(operator::SBPOperator) = operator.auxiliary_operators.Q
weighted_derivative_matrix(operator::FunctionSpaceSBPOperator) = operator.Q
gradient_matrix(operator::FunctionSpaceSBPOperator) = derivative_matrix(operator)

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
