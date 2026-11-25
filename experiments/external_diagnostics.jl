# Execute the upstream diagnostics function without running its full
# optimization and plotting campaign.
using LinearAlgebra, Printf
include("../test/fixtures/external_function_space.jl")

length(ARGS) == 1 || error("pass the path to the 2026_SBP_not_enough checkout")
source = read(joinpath(only(ARGS), "code", "operator_diagnostics.jl"), String)
start = findfirst("const NU = 1.0", source)
stop = findfirst("polybasis(d) =", source)
!isnothing(start) && !isnothing(stop) && first(start) < first(stop) ||
    error("upstream diagnostic definitions changed")
include_string(Main, source[first(start):prevind(source, first(stop))],
               "upstream operator_diagnostics.jl")

struct ReferenceOperator
    x::Vector{Float64}
    weights::Vector{Float64}
    D::Matrix{Float64}
end

grid(op::ReferenceOperator) = op.x
mass_matrix(op::ReferenceOperator) = Diagonal(op.weights)
Base.Matrix(op::ReferenceOperator) = copy(op.D)

for ref in EXTERNAL_FUNCTION_SPACE
    basis, derivatives = ref.name === :linear ? ((one,identity),(zero,one)) :
                         ((one,sinpi,cospi),(zero,x->π*cospi(x),x->-π*sinpi(x)))
    op = ReferenceOperator(ref.x, ref.weights, ref.D)
    report = diagnose(string(ref.name), op, basis, derivatives)
    @assert report.sv ≈ ref.singular_values
    @assert report.ranks[3] == (ref.name === :linear ? 3 : 6)
end
