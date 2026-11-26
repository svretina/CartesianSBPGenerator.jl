"""Construction-independent diagnostics for boundary-inclusive Cartesian SBP operators."""
module Diagnostics

using LinearAlgebra, SparseArrays
using ..CartesianSBPGenerator: derivative_matrix, mass_matrix, boundary_matrix,
    _rank_of_weighted_derivative

export OperatorData, diagnose_operator, algebraic_diagnostics, taylor_moments,
       norm_diagnostics, stencil_diagnostics, differentiation_errors,
       nullspace_diagnostics, spectral_diagnostics, sat_matrix, energy_transform,
       rk4_limit, rk4_polynomial, pseudospectrum, dispersion_diagnostics,
       bulk_stencil, variable_split

include("Algebraic.jl")
include("Spectrum.jl")
include("Dispersion.jl")

end
