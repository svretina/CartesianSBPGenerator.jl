@testset "analytical equal-order perturbation" begin
    for p in (2,4,6)
        baseline = design_sbp(2p,p; norm=:diagonal, skew_reconstruction=:analytical)
        base = instantiate(baseline; N=41, domain=(-1//2,3//2))
        op = CartesianSBPGenerator.Experimental.analytical_perturbation(
            baseline; N=41, domain=(-1//2,3//2),
        )
        @test op.H == base.H
        @test op.B == base.B
        K = (op.Q - base.auxiliary_operators.Q) / op.diagnostics.epsilon
        @test K == -K'
        @test all(iszero(K[i,j]) for i in 1:41, j in 1:41 if abs(i-j)>p+2)
        @test op.diagnostics.bandwidth_D <= max(baseline.closure_width-1,p+2)
        @test op.D == -reverse(op.D;dims=(1,2))
        @test op.H * op.D + op.D' * op.H == op.B
        @test op.diagnostics.nullity == 1
        X,Y = polynomial_matrices(op.x,p+1,1)
        @test all(iszero,K*X[:,1:(p+1)])
        @test op.D*X[:,1:(p+1)] == Y[:,1:(p+1)]
        @test all(!iszero,op.D*X[:,end]-Y[:,end])
        C = CartesianSBPGenerator.Experimental.leading_error_coefficients(op,p)
        i = op.diagnostics.closure_depth+1
        @test C[i] == -C[i+1]
        @test abs(C[i]) == 1//(12*2^(p-2))
        if p == 2
            @test op.diagnostics.bulk_taps[1] ==
                [1//192,0,7//96,-2//3,0,2//3,-7//96,0,-1//192]
            @test C[i] == -1//12
            @test_throws ArgumentError CartesianSBPGenerator.Experimental.analytical_perturbation(
                baseline; N=40,
            )
            @test_throws ArgumentError CartesianSBPGenerator.Experimental.analytical_perturbation(
                baseline; N=41, epsilon=0,
            )
        end
    end
end
