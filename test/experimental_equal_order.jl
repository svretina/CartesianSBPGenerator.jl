using CartesianSBPGenerator: Experimental

@testset "uniform polynomial element reference" begin
    for p in (2, 4, 6)
        op = Experimental.polynomial_element(p; domain=(-1//2, 3//2))
        D, H, B = derivative_matrix(op), mass_matrix(op), boundary_matrix(op)
        X, Y = polynomial_matrices(op.x, p + 1, 1)
        @test D * X[:, 1:(p + 1)] == Y[:, 1:(p + 1)]
        @test all(!iszero, D * X[:, end] - Y[:, end])
        @test H * D + transpose(D) * H == B
        @test all(k -> det(H[1:k, 1:k]) > 0, 1:(p + 1))
        @test op.diagnostics.nullity == 1

        element_error(family) = maximum(abs, family.D * family.x .^ (p + 1) - (p + 1) .* family.x .^ p)
        fine = Experimental.polynomial_element(p; domain=(-1//2, 1//2))
        @test element_error(op) == 2^p * element_error(fine)
    end
end

@testset "wider unknown stencil obstruction" begin
    for p in (2, 4, 6), extra_radius in (1, 2)
        result = Experimental.wider_stencil_obstruction(p; radius=p ÷ 2 + extra_radius,
                                                         closure_depth=p + extra_radius)
        @test result.status === :consistent
        @test result.forced_zero
        @test nullity(result.solution) > 0
    end
end

@testset "explicit periodic and compact experimental operators" begin
    cases = ((2, 2, 5, 2, 0, 17), (4, 2, 9, 4, 0, 25), (6, 2, 15, 6, 0, 37),
             (4, 1, 6, 3, 2, 17), (6, 1, 10, 4, 3, 25))
    for (p, period, depth, radius, bandwidth, n) in cases
        system = Experimental.assemble_stencil_experiment(p; N=n, period, closure_depth=depth,
                    stencil_radius=radius, norm_bandwidth=bandwidth)
        result = Experimental.search_stencil(system)
        @test result.status === :success
        op = result.operator
        refined = Experimental.refine_stencil(system, op; N=2n - 1, domain=(-1//2, 3//2))
        for candidate in (op, refined)
            D, H, Q = derivative_matrix(candidate), mass_matrix(candidate), weighted_derivative_matrix(candidate)
            X, Y = polynomial_matrices(candidate.x, p + 1, 1)
            @test D * X[:, 1:(p + 1)] == Y[:, 1:(p + 1)]
            @test all(!iszero, D * X[:, end] - Y[:, end])
            @test H * D == Q
            @test Q + transpose(Q) == boundary_matrix(candidate)
            @test H == reverse(H; dims=(1, 2))
            @test D == -reverse(D; dims=(1, 2))
            @test isposdef(Symmetric(Float64.(H)))
            @test candidate.diagnostics.nullity == 1
            @test candidate.diagnostics.bandwidth_H == bandwidth
        end

        hhat(candidate) = candidate.H / candidate.diagnostics.spacing
        @test hhat(op)[1:depth, 1:depth] == hhat(refined)[1:depth, 1:depth]
        @test op.Q[1:depth, 1:(depth + system.stencil_radius)] ==
              refined.Q[1:depth, 1:(depth + system.stencil_radius)]
        if period == 2
            data = Experimental.alternating_characteristics(system, op)
            @test data.weighted_skew
            @test data.parity_compatible
            @test data.weighted_moment == 0
            @test all(iszero, data.cut_flux)
            @test all(!iszero, data.leading)
            weighted_difference = data.weights[1] * data.moments[1][end] -
                                  data.weights[2] * data.moments[2][end]
            @test weighted_difference == (data.weights[1] - data.weights[2]) *
                  (-1)^(p ÷ 2 + 1) * prod(big.(2:2:p) .^ 2)
            # Check the radial bulk formula independently on both phases.
            r, h = 2//7, 1//70
            points = r .+ h .* data.offsets
            for phase in 1:2
                c = data.stencils[phase]
                for k in 0:(p-2)
                    @test sum(c .* points .^ (k+2)) / (h*r^2) == (k+2)*r^(k-1)
                end
                deficit = sum(c .* points .^ (p+1)) / (h*r^2) - (p+1)*r^(p-2)
                @test deficit == data.moments[phase][end] * h^p / r^2
            end
            W = Diagonal(Float64.(collect(data.weights)))
            for theta in (0.0, 0.1, 0.7, π/2)
                M = Experimental.bloch_symbol(data, theta)
                @test W*M + M'*W ≈ zeros(2,2) atol=2e-14
            end
            @test Experimental.bloch_symbol(data, 0.0) ≈ zeros(2,2) atol=2e-14
            @test_throws ArgumentError Experimental.refine_stencil(system, op; N=n+1)
            first_row = depth + 1
            stencil(i) = op.diagnostics.spacing .* op.D[i, (i - radius):(i + radius)]
            @test stencil(first_row) != stencil(first_row + 1)
            @test stencil(first_row) == stencil(first_row + 2)
        end
    end
end
