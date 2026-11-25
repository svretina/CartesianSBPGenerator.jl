@testset "generic centered stencils" begin
    @test minimum_centered_radius(1, 2) == 1
    @test minimum_centered_radius(2, 6) == 3
    @test minimum_centered_radius(6, 12) == 8
    @test_throws ArgumentError minimum_centered_radius(1, 3)
    @test_throws ArgumentError centered_stencil(
        derivative_order=4, accuracy_order=4, radius=2,
    )

    first_derivative = centered_stencil(
        derivative_order=1, accuracy_order=4,
    )
    @test coefficients(first_derivative) ==
          Rational{BigInt}[1//12, -2//3, 0, 2//3, -1//12]
    @test nullity(first_derivative) == 0

    second_derivative = centered_stencil(
        derivative_order=2, accuracy_order=6,
    )
    @test coefficients(second_derivative) ==
          Rational{BigInt}[2//180, -27//180, 270//180, -490//180,
                           270//180, -27//180, 2//180]

    # Exercise the complete range required by M01. Moment and parity
    # identities are checked independently of the construction routine.
    for derivative_order in 1:6, accuracy_order in 2:2:12
        stencil = centered_stencil(
            derivative_order=derivative_order,
            accuracy_order=accuracy_order,
        )
        c = coefficients(stencil)
        offsets = collect(stencil.offsets)
        for degree in 0:(derivative_order + accuracy_order - 1)
            expected = degree == derivative_order ? factorial(big(derivative_order)) : 0
            @test sum(c .* (offsets .^ degree)) == expected
        end
        @test all(
            c[stencil.radius + 1 - j] ==
            (-1)^derivative_order * c[stencil.radius + 1 + j]
            for j in 1:stencil.radius
        )
    end

    wider = centered_stencil(
        derivative_order=1, accuracy_order=4, radius=3,
    )
    @test nullity(wider) == 1
    @test length(coefficients(wider; parameters=[1//7])) == 7
    @test_throws DimensionMismatch coefficients(wider; parameters=[1, 2])
end

@testset "polynomial basis" begin
    X, Xdd = polynomial_matrices(collect(0:3), 3, 2)
    @test X[:, 4] == Rational{BigInt}[0, 1, 8, 27]
    @test Xdd[:, 1] == zeros(Rational{BigInt}, 4)
    @test Xdd[:, 3] == fill(Rational{BigInt}(2), 4)
    @test Xdd[:, 4] == Rational{BigInt}[0, 6, 12, 18]
end

@testset "sampled accuracy bases" begin
    nodes = 0:3
    X, Y = polynomial_matrices(nodes, 5, 1)
    full = sample_basis(PolynomialBasis(5), nodes)
    @test (full.X, full.Y) == (X, Y)
    sparse_basis = MonomialBasis([0, 1, 3, 5])
    selected = sample_basis(sparse_basis, nodes)
    @test selected.X == X[:, [1, 2, 4, 6]]
    @test selected.Y == Y[:, [1, 2, 4, 6]]
    @test selected.labels == ["x^0", "x^1", "x^3", "x^5"]
    @test eltype(evaluate_derivative(sparse_basis, nodes; arithmetic=Float64)) === Float64
    @test BasisSamples(sin.(nodes[:,:]), cos.(nodes[:,:]), ["sin"]).labels == ["sin"]
    @test_throws DimensionMismatch BasisSamples(X, Y[:, 1:2], full.labels)
end
