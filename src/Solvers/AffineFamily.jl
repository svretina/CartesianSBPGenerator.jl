"""Return coefficients for a member of an affine solution family."""
function coefficients(
    stencil::CenteredStencil{T};
    parameters::AbstractVector=Vector{T}(),
) where {T}
    nfree = nullity(stencil)
    if isempty(parameters)
        return copy(stencil.particular)
    end
    length(parameters) == nfree || throw(DimensionMismatch(
        "stencil has $nfree free parameters, got $(length(parameters))",
    ))
    return stencil.particular + stencil.nullspace * parameters
end

function evaluate(
    solution::AffineSolution{T},
    parameters::AbstractVector=Vector{T}(),
) where {T}
    isnothing(solution.particular) &&
        throw(ArgumentError("no certified affine solution to evaluate ($(solution.status))"))
    nfree = nullity(solution)
    if isempty(parameters)
        return copy(something(solution.particular))
    end
    length(parameters) == nfree || throw(DimensionMismatch(
        "solution has $nfree free parameters, got $(length(parameters))",
    ))
    return something(solution.particular) + solution.nullspace * parameters
end
