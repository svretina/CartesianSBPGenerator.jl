include("fixtures/external_classical.jl")
include("fixtures/external_function_space.jl")
include("fixtures/external_matlab_fsbp.jl")

@testset "upstream Mattsson–Nordström exact fixtures" begin
    for ref in EXTERNAL_CLASSICAL
        p, m, r = ref.order, ref.depth, ref.order ÷ 2
        spec = SBPSpec(interior_order=p,boundary_order=r,norm=:diagonal,closure_depth=m)
        direct = assemble_constraints(spec; method=:classical)
        reduced = assemble_constraints(spec; method=:reduced)
        z = closure_variables(direct, ref.H, ref.Q)
        @test direct.constraints.A * z == direct.constraints.b
        @test reduced.constraints.A * diag(ref.H) == reduced.constraints.b
        @test all(>(0), diag(ref.H))
        @test ref.Q[:, 1:m] + transpose(ref.Q[:, 1:m]) == Diagonal(vcat(-1, zeros(Int,m-1)))
        X, Y = polynomial_matrices(0:(m+r-1), r, 1)
        @test ref.Q * X == ref.H * Y[1:m, :]
        @test vec(ref.center) == coefficients(reduced.interior_stencil)
        family = design_sbp(p,r; norm=:diagonal)
        @test family.left_norm == ref.H
        if p <= 4
            @test family.left_Q == ref.Q
        else
            fixed_norm = reconstruct_skew(reduced, ref.H)
            target = z[(m+1):end]
            coordinates = solve_constraints(LinearConstraintSystem(
                fixed_norm.nullspace, target - fixed_norm.particular,
            ))
            @test evaluate(fixed_norm, evaluate(coordinates)) == target
        end
    end
end

@testset "upstream fixed-norm FSBP fixtures and diagnostics" begin
    for ref in EXTERNAL_FUNCTION_SPACE
        n = length(ref.x)
        H, D = Diagonal(ref.weights), ref.D
        B = Matrix(Diagonal(vcat(-1.0, zeros(n-2), 1.0)))
        Q = H * D
        @test all(>(0), ref.weights)
        @test D * ref.X ≈ ref.Y atol=1e-10
        @test Q + Q' ≈ B atol=1e-12

        # Match the externally supplied grid, basis and mass. Both
        # reconstructions may choose different admissible skew entries.
        target = H * ref.Y - B * ref.X / 2
        skew = CartesianSBPGenerator._analytical_skew(ref.X, target, :numerical)
        A = CartesianSBPGenerator._skew_block(evaluate(skew), n, Float64)
        Δ = Q - B/2 - A
        @test Δ + Δ' ≈ zeros(n,n) atol=1e-12
        @test Δ * ref.X ≈ zeros(size(ref.X)) atol=1e-10

        @test svdvals(D) ≈ ref.singular_values atol=1e-12
        nullity = n - CartesianSBPGenerator._rank_of_weighted_derivative(D)
        @test nullity == (ref.name === :linear ? 4 : 1)
        sat = copy(D)
        sat[1,1] += inv(ref.weights[1])
        actual = eigvals(sat)
        @test maximum(λ -> minimum(abs.(actual .- λ)), ref.sat_eigenvalues) < 1e-9
    end
end

@testset "original MATLAB FSBP references" begin
    for ref in EXTERNAL_MATLAB_FSBP
        n = length(ref.x)
        H = Diagonal(ref.weights)
        B = Matrix(Diagonal(vcat(-1.0, zeros(n-2), 1.0)))
        @test all(>(0), ref.weights)
        @test ref.D * ref.X ≈ ref.Y atol=1e-12
        @test H * ref.D ≈ ref.Q atol=1e-12
        @test ref.Q + ref.Q' ≈ B atol=1e-12
        skew = CartesianSBPGenerator._analytical_skew(ref.X, H * ref.Y - B * ref.X / 2, :numerical)
        A = CartesianSBPGenerator._skew_block(evaluate(skew), n, Float64)
        @test nullity(skew) == 0
        @test A + B/2 ≈ ref.Q atol=1e-12
        if ref.name === :poly
            # The quadratic column augments linear exactness. The unique
            # three-node solution must reproduce both matrices.
            op = generate_sbp(PolynomialBasis(2);N=n)
            @test Float64.(op.D) ≈ ref.D atol=1e-12
            @test Float64.(op.H) ≈ H atol=1e-12
        end
    end
end
