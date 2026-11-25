@testset "public construction methods" begin
    for (method, construction) in ((:classical, :direct), (:reduced, :reduced))
        op = generate_sbp(4, 2; N=17, method)
        @test op.diagnostics.method === method
        @test op.diagnostics.construction === construction
        @test weighted_derivative_matrix(op) == mass_matrix(op) * derivative_matrix(op)
    end
    classical = generate_sbp(4, 2; N=17, method=:classical, arithmetic=:numerical)
    @test !classical.diagnostics.exact
    @test classical.H * classical.D + transpose(classical.D) * classical.H ≈ classical.B
    @test_throws ArgumentError design_sbp(4, 2; method=:classical, construction=:reduced)
    @test_throws ArgumentError design_sbp(4, 2; method=:function_space)
end

@testset "fixed-grid function-space SBP" begin
    for (norm, m, bandwidth) in ((:diagonal, nothing, 1), (:restricted_full, 3, 6))
        op = generate_sbp(PolynomialBasis(1); N=7, domain=(-1//2, 5//2),
                          norm, closure_depth=m, bandwidth)
        D, H, B = derivative_matrix(op), mass_matrix(op), boundary_matrix(op)
        Q = weighted_derivative_matrix(op)
        @test op.diagnostics.method === :function_space
        @test op.diagnostics.nullity == 1
        @test gradient_matrix(op) === D
        @test D * op.x == ones(7)
        @test D * ones(Int, 7) == zeros(7)
        @test H * D == Q
        @test Q + transpose(Q) == B
        @test det(Q[1:6, 1:6]) != 0
        @test all(k -> det(H[1:k, 1:k]) > 0, 1:7)
        @test H == reverse(H; dims=(1, 2))
        @test D == -reverse(D; dims=(1, 2))
        @test all(iszero(Q[i, j]) for j in 1:7, i in 1:7 if abs(i - j) > bandwidth)
        @test all(iszero, H[1, 2:end])
        if norm === :diagonal
            @test isdiag(H)
        else
            @test H[4, 4] == 1//2
            @test all(iszero, H[1:3, 4:7])
        end
    end

    # The trigonometric basis is independent of the polynomial closure API.
    # Test exactness on the full grid, including rows far from either endpoint.
    x = range(0.0, 1.0; length=7)
    sine, cosine = sinpi.(2 .* x), cospi.(2 .* x)
    samples = BasisSamples(hcat(sine, cosine), 2π .* hcat(cosine, -sine), ["sin", "cos"])
    op = generate_sbp(samples; method=:function_space, arithmetic=:numerical, bandwidth=2)
    @test op.x == collect(x)
    @test op.D * samples.X ≈ samples.Y atol=1e-11
    @test op.H * op.D + transpose(op.D) * op.H ≈ op.B atol=1e-11
    @test rank(op.D) == length(x) - 1
    @test isposdef(Symmetric(op.H))
    @test all(iszero(op.Q[i, j]) for j in eachindex(x), i in eachindex(x) if abs(i - j) > 2)
    @test_throws DimensionMismatch generate_sbp(samples; N=8, arithmetic=:numerical)
    impossible = BasisSamples(zeros(Int, 5, 1), ones(Int, 5, 1), ["incompatible"])
    @test_throws ErrorException generate_sbp(impossible)
end
