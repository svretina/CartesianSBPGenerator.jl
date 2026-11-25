"""Construct an exact centered finite-difference stencil family.

The returned coefficients are dimensionless: an application on a grid with
spacing `h` must be divided by `h^derivative_order`. Wider-than-minimum
stencils expose their remaining degrees of freedom through `nullspace`.
"""
function centered_stencil(;
    derivative_order::Integer,
    accuracy_order::Integer,
    radius::Union{Nothing,Integer}=nothing,
    arithmetic::Type=BigRational,
)
    minimum = minimum_centered_radius(derivative_order, accuracy_order)
    selected_radius = isnothing(radius) ? minimum : radius
    system = moment_system(
        derivative_order, accuracy_order, selected_radius; arithmetic=arithmetic,
    )
    solution = solve_constraints(system; backend=:exact)
    solution.status === :inconsistent && error(
        "the centered-stencil constraints are inconsistent",
    )
    T = eltype(something(solution.particular))
    return CenteredStencil{T}(
        derivative_order,
        accuracy_order,
        selected_radius,
        -selected_radius:selected_radius,
        something(solution.particular),
        solution.nullspace,
    )
end
