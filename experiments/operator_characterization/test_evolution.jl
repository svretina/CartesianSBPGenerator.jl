using Test
include("characterize.jl")

@testset "manufactured evolution and SAT identities" begin
    family=design_sbp(4,2;skew_reconstruction=:analytical)
    op=OperatorData(instantiate(family;N=41);depth=family.closure_depth)
    taps=bulk_stencil(op)
    models=(advection_model(op), advection_model(op;ppw=8),
            wave_model(op), wave_model(op;absorbing=true),
            periodic_model(op,taps), variable_model(op;compatible=true),
            variable_model(op;compatible=false,degenerate=true))
    for model in models
        r=evolve(model,[0.0,0.125,0.25];max_steps=1000,accuracy_factor=0.025)
        @test r.explicit
        @test maximum(abs,r.power_residual) <= 1e-11
        @test r.temporal_error <= 1e-6
        @test r.exponential_defect <= 1e-10
        @test all(isfinite,r.states)
        @test all(>=(0),r.energy)
        @test size(r.states)==(model.n,3)
        for (u,t) in zip(eachcol(reduce(hcat,model.reference.(r.times))),r.times)
            @test model.power(u,model.M*u) ≈ model.expected_power(u,t) atol=3e-11
        end
    end
    @test periodic_model(op,taps).M ≈ -periodic_model(op,taps).M'
    @test !evolve(first(models),[0.0,0.25];max_steps=1).explicit
    family=design_sbp(6,4;norm=:restricted_full,skew_reconstruction=:analytical)
    num=OperatorData(instantiate(family;N=41);depth=family.closure_depth)
    a=1 .+ 0.5sinpi.(2Float64.(num.x))
    fixed=variable_split(num,a;compatible=true)
    standard=variable_split(num,a;compatible=false)
    @test fixed.commutator>0
    @test fixed.residual<1e-9
    @test standard.residual>1e-6
end

@testset "analytical template transfer" begin
    for p in (2,4,6)
        arguments=(;depth=max(design_sbp(2p,p;skew_reconstruction=:analytical).closure_depth,p+2),period=2)
        source=OperatorData(E.analytical_perturbation(p;N=41);arguments...)
        refined=transfer_local(source,65,(0//1,1//1))
        independent=OperatorData(E.analytical_perturbation(p;N=65);arguments...)
        @test refined.D==independent.D
        @test refined.H==independent.H
        @test_throws ArgumentError transfer_local(source,42,(0//1,1//1))
    end
end

@testset "local conventional refinement and checkpoint" begin
    for (i,p,norm) in ((4,2,:diagonal),(8,6,:restricted_full))
        family=design_sbp(i,p;norm,skew_reconstruction=:analytical)
        depth=family.closure_depth
        template=OperatorData(instantiate(family;N=2*family.closure_width+1);depth)
        domain=(-1//2,3//2)
        sparse_op=transfer_local(template,65,domain)
        reference=OperatorData(instantiate(family;N=65,domain);depth)
        @test sparse_op.D==reference.D
        @test sparse_op.H==reference.H
        @test sparse_op.B==reference.B
        large=transfer_local(template,1025,domain)
        @test nnz(large.D) <= 1025*(2*family.closure_width+1)
        @test algebraic_diagnostics(large;minimum_order=p,exact_rank_max=0).certified
    end
    mktempdir() do dir
        filename=joinpath(dir,"nested","checkpoint.jls")
        original=(;rows=[(n=41,error=0.1)],complete=true)
        checkpoint(filename,original)
        @test deserialize(filename)==original
        @test readdir(dirname(filename))==["checkpoint.jls"]
    end
end
