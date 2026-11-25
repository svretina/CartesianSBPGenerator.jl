include("certify_alternating.jl") # Common exact and PDE diagnostics.

function write_table(path, A)
    open(path, "w") do io
        for row in eachrow(A)
            println(io, join(row, '\t'))
        end
    end
end

function dump_matrices(dir, p, n, name, op)
    h = (last(op.x)-first(op.x))/(n-1)
    prefix = joinpath(dir, "p$(p)_n$(n)_$(name)")
    write_table(prefix * "_D.tsv", Float64.(h .* op.D))
    write_table(prefix * "_rows.tsv",
        hcat(Float64.(diag(op.H)/h), Float64.(Experimental.leading_error_coefficients(op,p))))
end

option = findfirst(arg -> startswith(arg,"--matrix-dir="), ARGS)
matrix_dir = isnothing(option) ? nothing : split(ARGS[option],"=";limit=2)[2]
isnothing(matrix_dir) || mkpath(matrix_dir)

for (p,depth,N) in DESIGNS
    system, original = seed(p,depth,N)
    baseline = design_sbp(2p,p; norm=:diagonal, skew_reconstruction=:analytical)
    for n in (41,61,81)
        direct = Experimental.refine_stencil(system,original;N=n)
        conventional = instantiate(baseline;N=n)
        analytical = Experimental.analytical_perturbation(baseline;N=n)
        h = 1//(n-1)
        families = (:direct,:analytical,:conventional)
        candidates = (direct,analytical,conventional)
        depths = (depth,analytical.diagnostics.closure_depth,baseline.closure_depth)
        sats = map(sat_matrices,candidates)
        common_step = min(0.01 * Float64(h),minimum(sat.max_step for sat in sats)/8)
        for (name,op,m,sat) in zip(families,candidates,depths,sats)
            C = Experimental.leading_error_coefficients(op,p)
            radius = maximum(abs(i-j) for i in 1:n,j in 1:n if !iszero(op.D[i,j]))
            nullity = n - CartesianSBPGenerator._rank_of_weighted_derivative(weighted_derivative_matrix(op))
            sat_radius = Float64(h)*maximum(abs,eigvals(sat.plus))
            println((; p,n,name,bandwidth=radius,nullity,sat_radius,
                     spectral_cfl=sat.max_step/Float64(h), nnz=count(!iszero,op.D),
                     leading_extrema=Float64.(extrema(C)), operator_metrics(op,m)...,
                     truncation=truncation(op,m)))
            if "--pde" in ARGS
                benchmark = pde_benchmark(sat;step=common_step)
                all(isfinite,(benchmark.advection_error,benchmark.wave_error,benchmark.afine,benchmark.wfine)) ||
                    error("nonfinite PDE comparison")
                println((; p,n,name,task=:pde,benchmark...))
            end
            isnothing(matrix_dir) || n != 81 || dump_matrices(matrix_dir,p,n,name,op)
        end
    end
end
