function closure_variables(system, H, Q)
    m = system.closure_depth
    entries = [H[i, j] for (i, j) in system.norm_entries]
    skew = [Q[i, j] for j in 2:m for i in 1:(j - 1)]
    return vcat(entries, skew)
end

@testset "reduced and direct closure systems" begin
    cases = ((2, 1, :diagonal), (4, 2, :diagonal), (6, 3, :diagonal),
             (8, 4, :diagonal), (4, 3, :restricted_full))
    for (p, b, norm) in cases
        spec = SBPSpec(interior_order=p, boundary_order=b, norm=norm)
        reduced = assemble_constraints(spec)
        direct = assemble_constraints(spec; construction=:direct)
        reduced_solution = solve_constraints(reduced)
        direct_solution = solve_constraints(direct)
        @test size(reduced.constraints.A, 2) == length(reduced.norm_entries)
        @test size(direct.constraints.A, 2) > size(reduced.constraints.A, 2) ||
              reduced.closure_depth == 1

        direct_family = design_sbp(p, b; norm, construction=:direct)
        for method in (:linear, :analytical)
            family = design_sbp(p, b; norm, skew_reconstruction=method)
            @test family.diagnostics.constraint_nullity == nullity(direct_solution)
            @test direct.constraints.A *
                  closure_variables(direct, family.left_norm, family.left_Q) == direct.constraints.b
            @test reduced.constraints.A *
                  [direct_family.left_norm[i, j] for (i, j) in reduced.norm_entries] ==
                  reduced.constraints.b

            skew = reconstruct_skew(reduced, family.left_norm; method)
            rank_X = solve_constraints(LinearConstraintSystem(
                reduced.X, zeros(eltype(reduced.X), reduced.closure_depth),
            )).rank
            codimension = reduced.closure_depth - rank_X
            @test nullity(skew) == codimension * (codimension - 1) ÷ 2
            @test nullity(reduced_solution) + nullity(skew) == nullity(direct_solution)

            if nullity(skew) > 0
                variables = closure_variables(direct, family.left_norm, family.left_Q)
                variables[(length(direct.norm_entries) + 1):end] .=
                    evaluate(skew, ones(Int, nullity(skew)))
                @test direct.constraints.A * variables == direct.constraints.b
            end
        end
    end

    # The symmetric compatibility equations alone admit the classical norm.
    # The extra derivative is orthogonal to the sampled polynomial columns
    # in that norm, but its zero value column requires H*Y_extra to vanish.
    spec = SBPSpec(interior_order=2, boundary_order=1, closure_depth=3)
    conflicting = BasisSamples(zeros(Int, 4, 1), reshape([2, -2, 1, 0], :, 1), ["conflict"])
    for construction in (:direct, :reduced)
        system = assemble_constraints(spec; construction, basis=conflicting)
        @test solve_constraints(system).status === :inconsistent
        if construction === :reduced
            classical_weights = [1//2, 1, 1]
            @test system.constraints.A[1:6, :] * classical_weights == system.constraints.b[1:6]
        end
    end
end

@testset "nonpolynomial closure exactness" begin
    m = 6
    values = exp.(-(0:m))
    samples = BasisSamples(reshape(values, :, 1), reshape(-values, :, 1), ["exp(-ξ)"])
    for method in (:linear, :analytical)
        family = design_sbp(2, 1; closure_depth=m, basis=samples, arithmetic=:numerical,
                            skew_reconstruction=method)
        @test all(>(0), family.left_weights)
        op = instantiate(family; N=25, domain=(0, 24))
        left, right = 1:m, (length(op.x) - m + 1):length(op.x)
        fleft = exp.(-op.x)
        fright = exp.(op.x .- op.x[end])
        @test (op.D * fleft)[left] ≈ -fleft[left] atol=1e-11
        @test (op.D * fright)[right] ≈ fright[right] atol=1e-11
        @test op.D * ones(length(op.x)) ≈ zeros(length(op.x)) atol=1e-11
        @test op.D * op.x ≈ ones(length(op.x)) atol=1e-11
        @test op.H * op.D + transpose(op.D) * op.H ≈ op.B atol=1e-11
    end
end
