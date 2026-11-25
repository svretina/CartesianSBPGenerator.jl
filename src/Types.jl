"""A linear equality system `A * z = b` with optional variable names."""
struct LinearConstraintSystem{T,MT<:AbstractMatrix{T},VT<:AbstractVector{T}}
    A::MT
    b::VT
    variable_names::Vector{Symbol}

    function LinearConstraintSystem(
        A::MT,
        b::VT,
        variable_names::AbstractVector{Symbol}=Symbol[],
    ) where {T,MT<:AbstractMatrix{T},VT<:AbstractVector{T}}
        size(A, 1) == length(b) || throw(DimensionMismatch(
            "A has $(size(A, 1)) rows but b has length $(length(b))",
        ))
        isempty(variable_names) || length(variable_names) == size(A, 2) ||
            throw(DimensionMismatch(
                "expected $(size(A, 2)) variable names, got $(length(variable_names))",
            ))
        names = isempty(variable_names) ?
            [Symbol(:z, i) for i in axes(A, 2)] : collect(variable_names)
        new{T,MT,VT}(A, b, names)
    end
end

"""The complete solution of a consistent linear system.

Every solution is `particular + nullspace * parameters`. For an inconsistent
system, `particular` is `nothing` and the nullspace has no columns.
"""
struct AffineSolution{T}
    status::Symbol
    rank::Int
    pivot_columns::Vector{Int}
    particular::Union{Nothing,Vector{T}}
    nullspace::Matrix{T}
end

nullity(solution::AffineSolution) = size(solution.nullspace, 2)

"""A centered finite-difference stencil and its possible free coefficients.

`particular` is one exact solution of the moment equations. Every stencil in
the family is obtained by adding `nullspace * parameters`.
"""
struct CenteredStencil{T}
    derivative_order::Int
    accuracy_order::Int
    radius::Int
    offsets::UnitRange{Int}
    particular::Vector{T}
    nullspace::Matrix{T}
end

nullity(stencil::CenteredStencil) = size(stencil.nullspace, 2)

function Base.show(io::IO, stencil::CenteredStencil)
    print(io, "CenteredStencil(d=$(stencil.derivative_order), ",
          "p=$(stencil.accuracy_order), radius=$(stencil.radius), ",
          "nullity=$(nullity(stencil)))")
end

abstract type AbstractSBPStrategy end
struct FirstDerivativeSBP <: AbstractSBPStrategy end
struct SecondDerivativeSBP <: AbstractSBPStrategy end

Base.@kwdef struct SBPSpec
    interior_order::Int
    boundary_order::Int
    derivative_order::Int = 1
    norm::Symbol = :auto
    order_policy::Symbol = :exact
    interior_radius::Union{Nothing,Int} = nothing
    closure_depth::Union{Nothing,Int} = nothing
    norm_bandwidth::Union{Nothing,Int} = nothing
    boundary_derivative_order::Union{Nothing,Int} = nothing
    compatibility::Symbol = :independent
    extra_accuracy::Symbol = :auto
    symmetric::Bool = true
end

"""A resolution-independent first-derivative SBP design."""
struct SBPFamily{T,ST<:CenteredStencil{T}}
    derivative_order::Int
    requested_interior_order::Int
    requested_boundary_order::Int
    achieved_interior_order::Int
    achieved_boundary_order::Int
    norm::Symbol
    closure_depth::Int
    closure_width::Int
    left_weights::Vector{T}
    left_norm::Matrix{T}
    left_Q::Matrix{T}
    interior_stencil::ST
    diagnostics::NamedTuple
end

"""An SBP operator instantiated on one uniform grid."""
struct SBPOperator{T,MT<:AbstractMatrix{T},VT<:AbstractVector{T},AT}
    derivative_order::Int
    requested_interior_order::Int
    requested_boundary_order::Int
    achieved_interior_order::Int
    achieved_boundary_order::Int
    x::VT
    D::MT
    H::MT
    B::MT
    auxiliary_operators::AT
    diagnostics::NamedTuple
end

"""A first-derivative SBP operator exact on a sampled space on one fixed grid."""
struct FunctionSpaceSBPOperator{T,DT<:NamedTuple}
    derivative_order::Int
    x::Vector{T}
    D::Matrix{T}
    H::Matrix{T}
    B::Matrix{T}
    Q::Matrix{T}
    basis_labels::Vector{String}
    diagnostics::DT
end

"""The outcome of one fixed closure configuration in the stencil search."""
struct SBPClosureAttempt
    norm::Symbol
    closure_depth::Int
    interior_radius::Int
    status::Symbol
end

"""A failed stencil construction, with the configurations that were tried."""
struct SBPConstructionError <: Exception
    interior_order::Int
    boundary_order::Int
    status::Symbol
    attempts::Vector{SBPClosureAttempt}
end

function Base.showerror(io::IO, err::SBPConstructionError)
    print(io, "SBP¹(", err.interior_order, ",", err.boundary_order,
          ") construction returned ", err.status)
    if err.status === :theoretical_obstruction
        print(io, ". Boundary exactness forces the interior moment of degree ",
              err.interior_order + 1, " to vanish, contradicting exact interior order ",
              err.interior_order, ". Use order_policy=:allow_upgrade to permit a higher interior order.")
    end
    if !isempty(err.attempts)
        last = err.attempts[end]
        print(io, " after ", length(err.attempts), " candidate(s); last was ", last.norm,
              " at depth ", last.closure_depth, " (", last.status, ")")
    end
end
