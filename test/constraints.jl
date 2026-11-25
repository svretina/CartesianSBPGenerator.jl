@testset "generic linear constraint engine" begin
    unique_system = LinearConstraintSystem([1 1; 1 -1], [4, 2], [:x, :y])
    exact = solve_constraints(unique_system)
    @test exact.status == :unique
    @test exact.rank == 2
    @test nullity(exact) == 0
    @test evaluate(exact) == Rational{BigInt}[3, 1]

    underdetermined = LinearConstraintSystem([1 2 3], [4])
    family = solve_constraints(underdetermined)
    @test family.status == :underdetermined
    @test family.rank == 1
    @test nullity(family) == 2
    z = evaluate(family, Rational{BigInt}[2, -1])
    @test underdetermined.A * z == underdetermined.b

    inconsistent = LinearConstraintSystem([1 1; 2 2], [1, 3])
    failure = solve_constraints(inconsistent)
    @test failure.status == :inconsistent
    @test isnothing(failure.particular)
    @test_throws ArgumentError evaluate(failure)

    numerical = solve_constraints(
        LinearConstraintSystem([1.0 1.0; 1.0 -1.0], [4.0, 2.0]);
        backend=:numerical,
    )
    @test numerical.status == :unique
    @test evaluate(numerical) ≈ [3.0, 1.0]
    @test_throws ArgumentError solve_constraints(unique_system; backend=:unknown)
end

@testset "numerical affine families" begin
    for (A, b) in (([1 2 3; 2 4 6], [1, 2]), ([1 2; 2 4; 3 6], [1, 2, 3]))
        solution = solve_constraints(LinearConstraintSystem(A, b); backend=:numerical)
        @test solution.status === :underdetermined
        @test solution.rank == 1
        @test nullity(solution) == size(A, 2) - 1
        @test A * evaluate(solution, ones(nullity(solution))) ≈ b
        @test norm(A * solution.nullspace) < 1e-12
    end
    inconsistent = LinearConstraintSystem([1.0 2.0; 2.0 4.0], [1.0, 3.0])
    @test solve_constraints(inconsistent; backend=:numerical).status === :inconsistent
    uncertain = LinearConstraintSystem([1.0 0; 0 1e-12], [1.0, 0])
    @test solve_constraints(uncertain; backend=:numerical).status === :numerically_uncertain
    empty = LinearConstraintSystem(zeros(0, 2), Float64[])
    @test nullity(solve_constraints(empty; backend=:numerical)) == 2
    @test nullity(solve_constraints(empty)) == 2
end

@testset "constraint input validation" begin
    @test_throws DimensionMismatch LinearConstraintSystem(zeros(2, 3), zeros(3))
    @test_throws DimensionMismatch LinearConstraintSystem(zeros(2, 3), zeros(2), [:x])
end
