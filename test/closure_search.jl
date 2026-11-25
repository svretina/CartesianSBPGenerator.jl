@testset "exact diagonal positivity feasibility" begin
    R = Rational{BigInt}
    A = R[1 0; 0 1; -1 -1]
    b = R[0, 0, 1]
    status, point = CartesianSBPGenerator._strict_feasibility(A, b, Ref(100))
    @test status === :feasible
    @test all(>(0), A * point + b)
    b[end] = 0
    @test first(CartesianSBPGenerator._strict_feasibility(A, b, Ref(100))) === :infeasible
    @test first(CartesianSBPGenerator._strict_feasibility(A, b, Ref(0))) === :unresolved
end

@testset "automatic closure search" begin
    for (p, b, depth) in ((10, 5, 11), (12, 6, 14))
        family = design_sbp(p, b; norm=:diagonal, interior_radius=:minimal, max_closure_depth=depth)
        @test family.closure_depth == depth
        @test (family.achieved_interior_order, family.achieved_boundary_order) == (p, b)
        @test all(>(0), family.left_weights)
        attempts = family.diagnostics.attempts
        @test [a.closure_depth for a in attempts] == collect((p ÷ 2):depth)
        @test attempts[end].status === :success
        @test attempts[end - 1].status === :certified_positivity_infeasible
        X, Y = polynomial_matrices(0:(family.closure_width - 1), b, 1)
        @test family.left_Q * X == family.left_norm * Y[1:depth, :]
    end

    restricted = design_sbp(6, 5; norm=:restricted_full, max_closure_depth=7)
    @test restricted.closure_depth == 7
    @test all(k -> det(restricted.left_norm[1:k, 1:k]) > 0, 1:7)
    @test all(iszero, restricted.left_norm[1, 2:end])
    op = instantiate(restricted; N=23)
    @test op.H * op.D + transpose(op.D) * op.H == op.B
    @test op.D * (op.x .^ 5) == 5 .* op.x .^ 4

    err = try
        design_sbp(4, 3; norm=:restricted_full, max_closure_depth=4)
    catch failure
        failure
    end
    @test err isa SBPConstructionError
    @test err.status === :search_limit_reached
    @test all(a -> a.status === :inconsistent_linear_system, err.attempts)
    @test [a.closure_depth for a in err.attempts] == [2, 3, 4]

    @test_throws SBPConstructionError design_sbp(2, 1; closure_depth=1, positivity_margin=1//2)
    @test design_sbp(2, 1; positivity_margin=1//3).left_weights == [1//2]
    @test_throws ArgumentError design_sbp(4, 2; max_closure_depth=1)
end

@testset "equal-order obstruction and explicit upgrade" begin
    for p in (2, 4, 6)
        err = try
            design_sbp(p, p; norm=:restricted_full)
        catch failure
            failure
        end
        @test err isa SBPConstructionError
        @test err.status === :theoretical_obstruction
        @test isempty(err.attempts)
        @test occursin("degree $(p + 1)", sprint(showerror, err))

        # Contract the exact reduced equations with the coefficient vector
        # of (y-x)^p, accounting for upper-triangular packing.
        for m in (p + 1, p + 2)
            spec = SBPSpec(interior_order=p, boundary_order=p,
                           norm=:restricted_full, closure_depth=m)
            system = assemble_constraints(spec)
            R = Rational{BigInt}
            witness = R[j + k == p ? (j == k ? 1 : 2) * (-1)^j * binomial(p, j) : 0
                        for k in 0:p for j in 0:k]
            stencil = system.interior_stencil
            moment = sum(coefficients(stencil) .* (big.(stencil.offsets) .^ (p + 1)))
            @test all(iszero, transpose(witness) * system.constraints.A)
            @test dot(witness, system.constraints.b) == moment
            @test !iszero(moment)
        end

        family = design_sbp(p, p; norm=:restricted_full, order_policy=:allow_upgrade,
                            max_closure_depth=12)
        @test (family.requested_interior_order, family.requested_boundary_order) == (p, p)
        @test (family.achieved_interior_order, family.achieved_boundary_order) == (p + 2, p)
        @test all(k -> det(family.left_norm[1:k, 1:k]) > 0, 1:family.closure_depth)
        op = instantiate(family; N=2 * family.closure_width + 1)
        @test op.H * op.D + transpose(op.D) * op.H == op.B
        @test op.D * op.x .^ p == p .* op.x .^ (p - 1)
    end
end
