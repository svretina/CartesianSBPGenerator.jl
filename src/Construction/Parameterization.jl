"""An affine matrix `constant + sum(parameters[i] * basis[i])`."""
struct MatrixParameterization{T}
    constant::Matrix{T}
    basis::Vector{Matrix{T}}
    parameter_names::Vector{Symbol}

    function MatrixParameterization(
        constant::Matrix{T},
        basis::Vector{Matrix{T}},
        parameter_names::AbstractVector{Symbol},
    ) where {T}
        length(basis) == length(parameter_names) || throw(DimensionMismatch(
            "basis and parameter_names must have the same length",
        ))
        all(B -> size(B) == size(constant), basis) ||
            throw(DimensionMismatch("all basis matrices must match the constant"))
        new{T}(constant, basis, collect(parameter_names))
    end
end

function evaluate(parameterization::MatrixParameterization{T}, parameters::AbstractVector) where {T}
    length(parameters) == length(parameterization.basis) ||
        throw(DimensionMismatch("incorrect number of matrix parameters"))
    result = copy(parameterization.constant)
    for (parameter, basis_matrix) in zip(parameters, parameterization.basis)
        result .+= parameter .* basis_matrix
    end
    return result
end

function dense_parameterization(rows::Integer, columns::Integer; arithmetic::Type=BigRational)
    rows > 0 && columns > 0 || throw(ArgumentError("matrix dimensions must be positive"))
    basis = Matrix{arithmetic}[]
    names = Symbol[]
    for j in 1:columns, i in 1:rows
        B = zeros(arithmetic, rows, columns)
        B[i, j] = one(arithmetic)
        push!(basis, B)
        push!(names, Symbol(:a_, i, :_, j))
    end
    return MatrixParameterization(zeros(arithmetic, rows, columns), basis, names)
end

function symmetric_parameterization(n::Integer; arithmetic::Type=BigRational)
    n > 0 || throw(ArgumentError("matrix dimension must be positive"))
    basis = Matrix{arithmetic}[]
    names = Symbol[]
    for j in 1:n, i in 1:j
        B = zeros(arithmetic, n, n)
        B[i, j] = one(arithmetic)
        B[j, i] = one(arithmetic)
        push!(basis, B)
        push!(names, Symbol(:s_, i, :_, j))
    end
    return MatrixParameterization(zeros(arithmetic, n, n), basis, names)
end

function skew_symmetric_parameterization(n::Integer; arithmetic::Type=BigRational)
    n > 0 || throw(ArgumentError("matrix dimension must be positive"))
    basis = Matrix{arithmetic}[]
    names = Symbol[]
    for j in 2:n, i in 1:(j - 1)
        B = zeros(arithmetic, n, n)
        B[i, j] = one(arithmetic)
        B[j, i] = -one(arithmetic)
        push!(basis, B)
        push!(names, Symbol(:a_, i, :_, j))
    end
    return MatrixParameterization(zeros(arithmetic, n, n), basis, names)
end
