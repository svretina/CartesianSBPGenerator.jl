@testset "restricted-full first derivative" begin
    family = design_sbp(4, 3; norm=:restricted_full)
    @test family.norm === :restricted_full
    @test family.closure_depth == 5
    @test family.diagnostics.constraint_nullity == 3
    @test (family.achieved_interior_order, family.achieved_boundary_order) == (4, 3)
    block = family.left_norm
    @test block == transpose(block)
    @test all(iszero, block[1, 2:end])
    @test !isdiag(block)
    @test all(k -> det(block[1:k, 1:k]) > 0, axes(block, 1))

    m = family.closure_depth
    for N in (2 * family.closure_width + 1, 29)
        op = instantiate(family; N, domain=(-1//2, 3//2))
        D, H, B = derivative_matrix(op), mass_matrix(op), boundary_matrix(op)
        left, interior, right = 1:m, (m + 1):(N - m), (N - m + 1):N
        h = op.diagnostics.spacing
        @test H[left, left] == h .* block
        @test H[interior, interior] == h .* Matrix(I, length(interior), length(interior))
        @test all(iszero, H[left, (m + 1):N])
        @test all(iszero, H[right, 1:(N - m)])
        @test H == reverse(H; dims=(1, 2))
        @test D == -reverse(D; dims=(1, 2))
        @test H * D == op.auxiliary_operators.Q
        @test H * D + transpose(D) * H == B
        for degree in 0:4
            expected = derivative_monomial.(op.x, degree, 1)
            residual = D * (op.x .^ degree) - expected
            @test all(iszero, residual[interior])
            if degree <= 3
                @test all(iszero, residual[left])
                @test all(iszero, residual[right])
            else
                @test any(!iszero, residual[left])
                @test any(!iszero, residual[right])
            end
        end
    end

    errors = map((41, 81, 161)) do N
        op = instantiate(family; N, domain=(0.0, 1.0))
        residual = abs.(op.D * exp.(op.x) - exp.(op.x))
        return [maximum(residual[rows]) for rows in (1:m, (m + 1):(N - m), (N - m + 1):N)]
    end
    for (coarse, fine) in zip(errors, Base.tail(errors))
        @test log2.(coarse ./ fine) ≈ [3, 4, 3] atol=0.15
    end

    @test design_sbp(4, 3).norm === :restricted_full
    @test design_sbp(4, 2; norm=:restricted_full, norm_bandwidth=0).left_norm ==
          design_sbp(4, 2; norm=:diagonal).left_norm
    @test_throws ArgumentError design_sbp(4, 3; norm=:restricted_full, norm_bandwidth=-1)
    shallow = SBPSpec(interior_order=4, boundary_order=3, norm=:restricted_full, closure_depth=4)
    @test solve_constraints(assemble_constraints(shallow)).status === :inconsistent
end
