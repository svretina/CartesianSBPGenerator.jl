"""Opt-in equal-order SBP research on fixed grids and polynomial elements."""
module Experimental

using LinearAlgebra
using ..CartesianSBPGenerator: BigRational, LinearConstraintSystem, AffineSolution,
    FunctionSpaceSBPOperator, SBPFamily, polynomial_matrices, solve_constraints, evaluate,
    design_sbp, instantiate, coefficients, weighted_derivative_matrix,
    _exact, _function_space_matrix, _positive_definite, _rank_of_weighted_derivative,
    _near_norm_solution

export polynomial_element, assemble_stencil_experiment, stencil_candidate,
       search_stencil, refine_stencil, wider_stencil_obstruction,
       alternating_characteristics, bloch_symbol, analytical_perturbation,
       leading_error_coefficients

include("EqualOrder.jl")
include("AlternatingAnalysis.jl")
include("Perturbation.jl")

end
