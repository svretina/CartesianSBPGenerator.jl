@testset "matrix parameterizations" begin
    dense = dense_parameterization(2, 3)
    @test length(dense.basis) == 6
    @test evaluate(dense, 1:6) == Rational{BigInt}[1 3 5; 2 4 6]

    symmetric = symmetric_parameterization(3)
    S = evaluate(symmetric, 1:6)
    @test S == transpose(S)
    @test length(symmetric.basis) == 6

    skew = skew_symmetric_parameterization(4)
    A = evaluate(skew, 1:6)
    @test A == -transpose(A)
    @test iszero(diag(A))
    @test length(skew.basis) == 6
end
