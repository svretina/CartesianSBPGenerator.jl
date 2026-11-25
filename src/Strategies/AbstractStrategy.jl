function _validate_spec(spec::SBPSpec)
    spec.interior_order > 0 || throw(ArgumentError("interior_order must be positive"))
    spec.boundary_order > 0 || throw(ArgumentError("boundary_order must be positive"))
    spec.derivative_order > 0 || throw(ArgumentError("derivative_order must be positive"))
    iseven(spec.interior_order) || throw(ArgumentError(
        "the initial centered-stencil backend requires an even interior_order",
    ))
    spec.norm in (:auto, :diagonal, :restricted_full) || throw(ArgumentError(
        "norm must be :auto, :diagonal, or :restricted_full",
    ))
    spec.order_policy in (:exact, :allow_upgrade) || throw(ArgumentError(
        "order_policy must be :exact or :allow_upgrade",
    ))
    spec.compatibility in (:independent, :compatible, :fully_compatible) ||
        throw(ArgumentError("unsupported compatibility mode $(spec.compatibility)"))
    return spec
end

function _strategy(spec::SBPSpec)
    spec.derivative_order == 1 && return FirstDerivativeSBP()
    spec.derivative_order == 2 && return SecondDerivativeSBP()
    throw(ArgumentError(
        "no validated SBP construction backend exists for derivative order " *
        "$(spec.derivative_order); centered_stencil supports it independently",
    ))
end
