@testset "construction-independent SBP diagnostics" begin
    using CartesianSBPGenerator.Diagnostics
    H = Rational{BigInt}[1//4 0 0; 0 1//2 0; 0 0 1//4]
    D = Rational{BigInt}[-2 2 0; -1 0 1; 0 -2 2]
    B = Matrix(Diagonal([-1,0,1]))
    op = OperatorData(D,H,B,[0//1,1//2,1//1];depth=1)
    result = diagnose_operator(op;degree=6)
    @test result.algebra.certified
    @test result.algebra.exact_rank == 2
    @test result.algebra.sbp_relative == 0
    @test result.algebra.taylor.orders == [1,2,1]
    @test result.algebra.taylor.leading == [1//2,1//6,-1//2]
    @test result.norms.eigenvalues == [0.5,0.5,1]
    @test result.norms.quadrature[1:3] == [0,0,1//24]
    @test result.norms.bilinear[2,2] == 1//24
    @test result.stencil.D == (;lower=1,upper=1)
    @test result.stencil.counts == [2,2,2]
    @test result.stencil.signed_displacement == [0.5,0,-0.5]
    @test result.kernel.numerical_rank == 2
    @test result.spectral.strict == :strict
    @test result.spectral.energy_abscissa <= 1e-14
    @test all(<=(1+1e-13),result.spectral.transient)
    @test !diagnose_operator(OperatorData(2D,H,B,op.x)).algebra.certified
    @test !diagnose_operator(OperatorData(D,-H,B,op.x)).algebra.certified
    @test_throws ArgumentError OperatorData(D,H,B,[0.0,0.4,1.0])
    numerical = OperatorData(Float64.(D),Float64.(H),Float64.(B),Float64.(op.x);depth=1)
    @test algebraic_diagnostics(numerical).certified
    @test rk4_limit(im.*[1.0,-1.0]) ≈ 2sqrt(2)
    @test rk4_limit([1.0,-1.0]) == 0
    @test abs(rk4_polynomial(-rk4_limit([-2.0])*2)) ≈ 1
    ps = pseudospectrum(Diagonal([-1.0,-2.0]),[-1.0,0.0],[0.0,1.0])
    @test ps ≈ [0 1; 1 sqrt(2)]
    split = variable_split(op,[0.75,1.0,1.25])
    @test split.commutator == 0
    @test split.residual <= 1e-14
    @test split.abscissa <= 1e-14
end

@testset "Fourier and Bloch diagnostics" begin
    using CartesianSBPGenerator.Diagnostics
    family = design_sbp(2,1)
    op = OperatorData(instantiate(family;N=17);depth=family.closure_depth)
    d = dispersion_diagnostics(op;samples=65)
    @test vec(d.frequencies) ≈ sin.(d.thetas)
    @test vec(d.group) ≈ cos.(d.thetas)
    @test vec(d.branch_energy) ≈ ones(65)
    for p in (2,4,6)
        baseline = design_sbp(2p,p;skew_reconstruction=:analytical)
        candidate = Experimental.analytical_perturbation(baseline;N=41)
        op = OperatorData(candidate;depth=candidate.diagnostics.closure_depth,period=2)
        d = dispersion_diagnostics(op;samples=33)
        e = Float64(candidate.diagnostics.epsilon)
        k(t)=sum(Float64(c)*sin(j*t) for (j,c) in
                 zip(baseline.interior_stencil.offsets,coefficients(baseline.interior_stencil)))
        for i in eachindex(d.thetas)
            t = d.thetas[i]
            beta = (-1)^(p÷2+1)*2.0^(p+2)*cos(t)*sin(t)^(p+1)
            analytical = eigvals(Symmetric([k(t) e*beta; e*beta k(t+π)]))
            @test sort(d.frequencies[i,:]) ≈ analytical atol=2e-14
        end
        @test vec(sum(d.branch_energy;dims=2)) ≈ ones(33)
        @test d.skew_residual <= 2e-14
    end
end
