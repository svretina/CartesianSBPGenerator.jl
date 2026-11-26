# Extract the norm and weighted derivative instead of truncating a dense H⁻¹Q.
function bulk_stencil(op::OperatorData)
    H = op.H ./ op.h
    Q = dropzeros!(H*(op.h .* op.D))
    bulk = row_groups(op).bulk
    function taps(M,i)
        ids, values = findnz(M[i,:])
        return Dict(j-i => v for (j,v) in zip(ids,values))
    end
    same(a,b) = a == b || (eltype(op.D) <: AbstractFloat &&
        all(k -> isapprox(get(a,k,0),get(b,k,0); atol=256eps(eltype(op.D)),rtol=256eps(eltype(op.D)),
            ), union(keys(a),keys(b))))
    return map(1:op.period) do phase
        rows = filter(i -> mod1(i,op.period)==phase,bulk)
        length(rows) >= 2 || throw(ArgumentError("need two bulk rows per phase"))
        norm, derivative = taps(H,first(rows)), taps(Q,first(rows))
        all(i -> same(taps(H,i),norm) && same(taps(Q,i),derivative),rows) ||
            throw(ArgumentError("interior does not repeat at the requested closure depth and period"))
        (; norm, derivative)
    end
end

function symbol(taps, theta::T, derivative=false) where {T<:AbstractFloat}
    period = length(taps)
    H, Q = zeros(Complex{T},period,period),zeros(Complex{T},period,period)
    for r in 1:period, (A, coefficients) in ((H,taps[r].norm),(Q,taps[r].derivative))
        for k in sort!(collect(keys(coefficients)))
            term = T(coefficients[k])*cis(k*theta)
            A[r,mod1(r+k,period)] += derivative ? im*k*term : term
        end
    end
    return H,Q
end

function dispersion_diagnostics(taps::AbstractVector; samples=4097,
                               thresholds=(0.1,0.01,0.001), min_theta=1e-4)
    period = length(taps)
    period in (1,2) && samples >= 3 || throw(ArgumentError("expected period 1 or 2 and at least three samples"))
    0 < min_theta < π/period || throw(ArgumentError("min_theta lies outside the positive Brillouin zone"))
    thetas = collect(range(min_theta,π/period;length=samples))
    frequencies, group = (zeros(samples,period) for _ in 1:2)
    branch_energy = zeros(samples,period)
    polarization = zeros(samples,period)
    norm_min = Inf
    skew_residual = 0.0
    previous = nothing
    for (i,theta) in enumerate(thetas)
        H,Q = symbol(taps,theta)
        Hp,Qp = symbol(taps,theta,true)
        residual=max(norm(Q+Q',Inf),norm(H-H',Inf))
        residual <= 512eps(Float64)*max(1,norm(H,Inf),norm(Q,Inf)) ||
            throw(ArgumentError("bulk symbols do not satisfy weighted skew symmetry"))
        norm_min = min(norm_min,eigmin(Hermitian(H)))
        skew_residual = max(skew_residual,residual)
        U = cholesky(Hermitian(H)).U
        eig = eigen(Hermitian(-im*(U'\Q/U)))
        vectors = eig.vectors
        smooth = U*ones(period)
        fractions = abs2.(vectors'*smooth)./sum(abs2,smooth)
        ids = if period == 1
            (1,)
        elseif isnothing(previous)
            first_branch = argmax(fractions)
            (first_branch,3-first_branch)
        elseif abs2(dot(previous[:,1],vectors[:,1])) + abs2(dot(previous[:,2],vectors[:,2])) >=
               abs2(dot(previous[:,1],vectors[:,2])) + abs2(dot(previous[:,2],vectors[:,1]))
            (1,2)
        else
            (2,1)
        end
        for (j,idx) in enumerate(ids)
            k = eig.values[idx]
            v = U\vectors[:,idx]
            frequencies[i,j] = k
            group[i,j] = real(dot(v,(-im*Qp-k*Hp)*v))
            branch_energy[i,j] = fractions[idx]
            polarization[i,j] = abs2(v[1])/sum(abs2,v)
        end
        previous = vectors[:,collect(ids)]
    end
    phase = frequencies ./ thetas
    function ppw(errors,tol)
        bad = findfirst(>(tol),errors)
        return isnothing(bad) ? 2π/last(thetas) : bad==1 ? Inf : 2π/thetas[bad-1]
    end
    crossings = [count(i -> frequencies[i,j]*frequencies[i+1,j] < 0,1:samples-1) for j in 1:period]
    # Report resolved sign-changing roots and separately the edge and origin.
    H0,Q0 = symbol(taps,0.0)
    U0 = cholesky(Hermitian(H0)).U
    origin = eigvals(Hermitian(-im*(U0'\Q0/U0)))
    return (; thetas, frequencies,phase,group,branch_energy,polarization,norm_min,
             skew_residual,crossings,origin,thresholds,ppw_phase=map(t -> ppw(abs.(phase[:,1].-1),t),thresholds),
             ppw_group=map(t -> ppw(abs.(group[:,1].-1),t),thresholds))
end

dispersion_diagnostics(op::OperatorData; kwargs...) =
    dispersion_diagnostics(bulk_stencil(op);kwargs...)
