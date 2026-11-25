@testset "diagonal-norm first derivatives" begin
    established_weights = Dict(
        (2, 1) => Rational{BigInt}[1//2],
        (4, 2) => Rational{BigInt}[17//48, 59//48, 43//48, 49//48],
        (6, 3) => Rational{BigInt}[
            13649//43200, 12013//8640, 2711//4320,
            5359//4320, 7877//8640, 43801//43200,
        ],
        (8, 4) => Rational{BigInt}[
            1498139//5080320, 1107307//725760, 20761//80640,
            1304999//725760, 299527//725760, 103097//80640,
            670091//725760, 5127739//5080320,
        ],
    )

    for ((interior_order, boundary_order), expected_weights) in established_weights
        family = design_sbp(interior_order, boundary_order; norm=:diagonal)
        @test family.left_weights == expected_weights
        @test all(>(0), family.left_weights)
        @test family.requested_interior_order == interior_order
        @test family.requested_boundary_order == boundary_order
        @test family.achieved_interior_order == interior_order
        @test family.achieved_boundary_order == boundary_order

        N = max(41, 2 * family.closure_width + 1)
        operator = instantiate(family; N=N, domain=(0, 2))
        D, H, B = derivative_matrix(operator), mass_matrix(operator), boundary_matrix(operator)
        Q = operator.auxiliary_operators.Q
        @test gradient_matrix(operator) === D
        @test H * D == Q
        @test Q + transpose(Q) == B
        @test isposdef(Symmetric(H))

        x = operator.x
        m = family.closure_depth
        for degree in 0:boundary_order
            exact = [derivative_monomial(node, degree, 1) for node in x]
            approximation = D * (x .^ degree)
            @test approximation[1:m] == exact[1:m]
            @test approximation[(N - m + 1):N] == exact[(N - m + 1):N]
        end
        for degree in 0:interior_order
            exact = [derivative_monomial(node, degree, 1) for node in x]
            approximation = D * (x .^ degree)
            @test approximation[(m + 1):(N - m)] == exact[(m + 1):(N - m)]
        end
    end
end

@testset "first-derivative request validation" begin
    @test_throws ArgumentError design_sbp(4, 3; norm=:diagonal)
    @test_throws ArgumentError design_sbp(4, 2; derivative_order=2)

    family = design_sbp(4, 2)
    @test_throws ArgumentError instantiate(family; N=5)
    operator = generate_sbp(2, 1; N=11)
    @test_throws ArgumentError laplacian_matrix(operator)

    floating = generate_sbp(2, 1; N=11, domain=(0.0, 1.0))
    @test eltype(floating.D) == Float64
    @test floating.H * floating.D + transpose(floating.D) * floating.H ≈ floating.B
end
