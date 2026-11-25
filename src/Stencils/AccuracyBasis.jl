"""Sampled basis values and first derivatives, with one function per column."""
struct BasisSamples{MX<:AbstractMatrix,MY<:AbstractMatrix}
    X::MX
    Y::MY
    labels::Vector{String}

    function BasisSamples(X::MX, Y::MY, labels) where {MX<:AbstractMatrix,MY<:AbstractMatrix}
        size(X) == size(Y) || throw(DimensionMismatch("basis values and derivatives must match"))
        length(labels) == size(X, 2) ||
            throw(DimensionMismatch("expected one label per basis column"))
        return new{MX,MY}(X, Y, collect(String, labels))
    end
end

"""The monomials of degrees zero through `degree`."""
struct PolynomialBasis
    degree::Int
    function PolynomialBasis(degree::Integer)
        degree >= 0 || throw(ArgumentError("degree must be nonnegative"))
        return new(degree)
    end
end

"""Monomials in the caller's chosen order."""
struct MonomialBasis
    degrees::Vector{Int}
    function MonomialBasis(degrees::AbstractVector{<:Integer})
        !isempty(degrees) && all(>=(0), degrees) && allunique(degrees) ||
            throw(ArgumentError("monomial degrees must be nonnegative, distinct, and nonempty"))
        return new(collect(Int, degrees))
    end
end

_degrees(basis::PolynomialBasis) = 0:basis.degree
_degrees(basis::MonomialBasis) = basis.degrees

function _evaluate_basis(basis, nodes, derivative_order, ::Type{T}) where {T}
    return T[derivative_monomial(node, degree, derivative_order; arithmetic=T)
             for node in nodes, degree in _degrees(basis)]
end

function evaluate(basis::Union{PolynomialBasis,MonomialBasis}, nodes::AbstractVector;
                  arithmetic::Type=BigRational)
    return _evaluate_basis(basis, nodes, 0, arithmetic)
end

function evaluate_derivative(basis::Union{PolynomialBasis,MonomialBasis}, nodes::AbstractVector;
                             arithmetic::Type=BigRational)
    return _evaluate_basis(basis, nodes, 1, arithmetic)
end

function sample_basis(basis::Union{PolynomialBasis,MonomialBasis}, nodes::AbstractVector;
                      arithmetic::Type=BigRational)
    return BasisSamples(evaluate(basis, nodes; arithmetic),
                        evaluate_derivative(basis, nodes; arithmetic),
                        ["x^$degree" for degree in _degrees(basis)])
end
