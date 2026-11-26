function energy_transform(H, L)
    U = cholesky(Symmetric(Matrix{Float64}(H))).U
    return U * Matrix{Float64}(L) / U
end

function sat_matrix(op::OperatorData; penalty=1.0, speed=1.0)
    speed > 0 && penalty >= 0 || throw(ArgumentError("nonnegative penalty and positive speed required"))
    H, D = Matrix{Float64}(op.H), Matrix{Float64}(op.D)
    e = zeros(length(op.x))
    e[1] = 1
    return -speed .* (D + penalty .* (Symmetric(H) \ e) * e')
end

rk4_polynomial(z) = 1 + z*(1 + z*(1/2 + z*(1/6 + z/24)))

function rk4_limit(values::AbstractVector)
    rho = maximum(abs,values)
    iszero(rho) && return Inf
    tol = 128eps(Float64)*rho
    maximum(real,values) > tol && return 0.0
    bound = Inf
    for lambda in values
        iszero(lambda) && continue
        ray = complex(min(real(lambda),0),imag(lambda))/abs(lambda)
        lo, hi = 0.0, 4.0
        # Find the first exit along the ray before refining the crossing.
        for t in range(1/16,4;length=64)
            if abs2(rk4_polynomial(t*ray)) > 1+8eps(Float64)
                hi = t
                break
            end
            lo = t
        end
        for _ in 1:44
            mid = (lo+hi)/2
            if abs2(rk4_polynomial(mid*ray)) <= 1+8eps(Float64)
                lo = mid
            else
                hi = mid
            end
        end
        bound = min(bound,lo/abs(lambda))
    end
    return bound
end

function rk4_matrix(M, step)
    Z = step .* M
    Z2 = Z*Z
    return I + Z + Z2*(I/2 + Z/6 + Z2/24)
end

function nullspace_diagnostics(op::OperatorData; rtol=nothing)
    D = Matrix{Float64}(op.h .* op.D)
    H = Matrix{Float64}(op.H)
    F = svd(D)
    U = cholesky(Symmetric(H)).U
    S = svdvals(U * D / U)
    n = length(op.x)
    relative = something(rtol,n*eps(Float64))
    numerical_rank = count(>(relative*first(F.S)),F.S)
    noise = (-1.0).^(0:n-1)
    window = sinpi.(Float64.((op.x.-first(op.x))./(last(op.x)-first(op.x)))).^4
    ratio(v,y) = norm(U*y)/norm(U*v)
    combined = hcat(D*noise,D*(noise.*window)) ./ Float64(op.h)
    return (; singular=F.S, energy_singular=S, numerical_rank,
        ranks=[count(>(t*first(F.S)),F.S) for t in (1e-6,1e-8,1e-10,1e-12,1e-14,1e-16)],
        sigma_nonzero=F.S[n-1], physical_sigma_nonzero=F.S[n-1]/Float64(op.h),
        condition_nonzero=F.S[1]/F.S[n-1],
        near_vectors=F.V[:,max(1,n-2):n],
        checkerboard=ratio(noise,combined[:,1]),
        windowed_checkerboard=ratio(noise.*window,combined[:,2]), checkerboard_derivative=combined)
end

function spectral_diagnostics(op::OperatorData; sat_penalty=1.0, safety=0.7,
                              steps=(1,8,32), times=(0.0,0.25,0.5,1.0,2.0),
                              advanced=true)
    h = Float64(op.h)
    raw = eigvals(Matrix{Float64}(op.h .* op.D))
    L = h .* sat_matrix(op;penalty=sat_penalty)
    F = eigen(L)
    A = energy_transform(op.H,L)
    U = cholesky(Symmetric(Matrix{Float64}(op.H))).U
    modes = U * F.vectors
    m = op.depth
    edge = vcat(1:m,length(op.x)-m+1:length(op.x))
    localization = vec(sum(abs2,modes[edge,:];dims=1)./sum(abs2,modes;dims=1))
    radius, abscissa = maximum(abs,F.values),maximum(real,F.values)
    tol = 128eps(Float64)*radius
    strict = abscissa < -tol ? :strict : abscissa > tol ? :unstable : :unresolved
    cfl = rk4_limit(F.values)
    abscissa_energy = eigmax(Symmetric((A+A')/2))
    departure = norm(A'*A-A*A')/max(norm(A)^2,eps(Float64))
    if advanced && isfinite(cfl)
        step = safety*cfl
        G = rk4_matrix(A,step)
        amplification = [opnorm(G^k,2) for k in steps]
        transient = [opnorm(exp(t*A),2) for t in times]
        # A sampled contractivity interval with a roundoff allowance.
        trial = range(0,cfl;length=33)
        failed = findfirst(v -> opnorm(rk4_matrix(A,v),2) > 1+1e-10,trial)
        energy_cfl = isnothing(failed) ? cfl : trial[max(1,failed-1)]
        numerical_range = numerical_range_boundary(A)
    else
        amplification, transient, energy_cfl = Float64[],Float64[],missing
        numerical_range = ComplexF64[]
    end
    return (; raw, raw_radius=maximum(abs,raw),raw_abscissa=maximum(real,raw),
         sat=F.values, radius, abscissa, strict, tol, localization,
         most_localized=F.vectors[:,argmax(localization)], cfl, sat_penalty,
         energy_abscissa=abscissa_energy, departure, steps, safety, amplification,
         times, transient, energy_cfl, numerical_range, energy_matrix=A)
end

function numerical_range_boundary(A; samples=49)
    return map(range(0,2π;length=samples)) do angle
        rotated=cis(-angle).*A
        projection=Hermitian((rotated+rotated')/2)
        v=only(eachcol(eigen(projection,size(A,1):size(A,1)).vectors))
        dot(v,A*v)
    end
end

function pseudospectrum(A, real_axis, imag_axis)
    return [last(svdvals(complex(z,y)*I-A)) for z in real_axis,y in imag_axis]
end

function variable_split(op::OperatorData, coefficient; compatible=true, penalty=1.0)
    a = Float64.(coefficient)
    all(>=(0),a) || throw(ArgumentError("this inflow treatment requires nonnegative speeds"))
    H = Matrix{Float64}(op.H)
    D = Matrix{Float64}(op.D)
    A = Diagonal(a)
    Q = H*D
    e = zeros(length(a))
    e[1] = 1
    lift = Symmetric(H)\e
    L = compatible ? -0.5 .* (Symmetric(H)\(A*Q+Q*A)) : -0.5 .* (A*D+D*A)
    L -= penalty*a[1] .* (lift*e')
    flux = -Matrix{Float64}(op.B)*A - 2penalty*a[1] .* (e*e')
    transformed = energy_transform(H,L)
    return (; L, lift, penalty, commutator=opnorm(H*A-A*H,2)/
            max(opnorm(H,2)*maximum(a),eps(Float64)),
            residual=opnorm(H*L+L'*H-flux,Inf),
            abscissa=eigmax(Symmetric((transformed+transformed')/2)))
end
