"""Uniform-grid derivative, norm, and boundary matrices. Matrices use physical units."""
struct OperatorData{T}
    D::SparseMatrixCSC{T,Int}
    H::SparseMatrixCSC{T,Int}
    B::SparseMatrixCSC{T,Int}
    x::Vector{T}
    h::T
    depth::Int
    period::Int
end

function OperatorData(D::AbstractMatrix, H::AbstractMatrix, B::AbstractMatrix, x;
                      depth::Integer=0, period::Integer=1)
    n = length(x)
    n >= 2 && size(D) == size(H) == size(B) == (n,n) ||
        throw(DimensionMismatch("D, H, B must match a grid with at least two nodes"))
    0 <= 2depth < n && period in (1,2) ||
        throw(ArgumentError("invalid closure depth or bulk period"))
    T = promote_type(eltype(D), eltype(H), eltype(B), eltype(x))
    T = T <: Union{Integer,Rational} ? Rational{BigInt} : T
    T <: Real || throw(ArgumentError("real matrices and nodes are required"))
    nodes = T.(collect(x))
    h = (last(nodes) - first(nodes))/(n-1)
    all(isfinite,nodes) && h > 0 || throw(ArgumentError("invalid grid"))
    reference = first(nodes) .+ (0:n-1).*h
    uniform = T <: Rational ? nodes == reference :
        all(abs.(nodes-reference) .<= 16eps(T)*max(one(T),maximum(abs,nodes)))
    uniform || throw(ArgumentError("expected a uniform grid"))
    matrices = map(M -> dropzeros!(copy(SparseMatrixCSC{T,Int}(sparse(M)))), (D,H,B))
    all(M -> all(isfinite,nonzeros(M)),matrices) || throw(ArgumentError("nonfinite matrix"))
    return OperatorData(matrices...,nodes,h,Int(depth),Int(period))
end

OperatorData(op; kwargs...) =
    OperatorData(derivative_matrix(op),mass_matrix(op),boundary_matrix(op),op.x; kwargs...)

function bandwidth(M)
    I,J,_ = findnz(M)
    return (lower=maximum(I-J;init=0), upper=maximum(J-I;init=0))
end

function exact_positive(H)
    width = bandwidth(H).lower
    iszero(width) && return all(>(0),diag(H))
    A = Matrix(H)
    for i in axes(A,1)
        pivot = A[i,i]
        pivot > 0 || return false
        rows = i+1:min(i+width,size(A,1))
        for j in rows, k in rows
            A[j,k] -= A[j,i] * A[i,k] / pivot
        end
    end
    return true
end

function taylor_moments(Dhat; degree::Integer=28, rtol=nothing)
    degree >= 1 || throw(ArgumentError("degree must be positive"))
    T = eltype(Dhat)
    exact = T <: Rational
    tol = isnothing(rtol) ? (exact ? 0 : 64eps(T)) : rtol
    n = size(Dhat,1)
    C = zeros(T,n,degree+1)
    scale = zeros(T,n,degree+1)
    I,J,V = findnz(sparse(Dhat))
    for (i,j,value) in zip(I,J,V)
        term = value
        for k in 0:degree
            C[i,k+1] += term
            scale[i,k+1] += abs(term)
            term *= T(j-i)/T(k+1)
        end
    end
    C[:,2] .-= 1
    scale[:,2] .+= 1
    orders = map(axes(C,1)) do i
        first_bad = findfirst(k -> abs(C[i,k]) > tol*scale[i,k], axes(C,2))
        isnothing(first_bad) ? degree : first_bad-2
    end
    leading = [p < degree ? C[i,p+2] : missing for (i,p) in enumerate(orders)]
    return (; coefficients=C, orders, leading, degree, rtol=tol)
end

function algebraic_diagnostics(op::OperatorData; degree=28, minimum_order=1,
                               exact_order=false, reflection=true, exact_rank_max=129)
    n = length(op.x)
    Dhat,Hhat = op.h .* op.D, op.H ./ op.h
    Q = dropzeros!(Hhat*Dhat)
    E = Q + Q' - op.B
    DR = Dhat + reverse(Dhat;dims=(1,2))
    HR = Hhat - reverse(Hhat;dims=(1,2))
    exact = eltype(op.D) <: Rational
    tol = exact ? 0 : 128eps(eltype(op.D))
    scaled_residual(A,scale) = opnorm(A,Inf)/max(scale,one(scale))
    sbp_absolute = opnorm(E,Inf)
    scale = 2opnorm(Q,Inf)+opnorm(op.B,Inf)
    sbp_relative = iszero(scale) ? zero(sbp_absolute) : sbp_absolute/scale
    taylor = taylor_moments(Dhat;degree)
    symmetric = scaled_residual(Hhat-Hhat',opnorm(Hhat,Inf)) <= tol
    positive = symmetric && (exact ? exact_positive(Hhat) :
                                              isposdef(Symmetric(Matrix{Float64}(Hhat))))
    exact_rank = exact && n <= exact_rank_max ? _rank_of_weighted_derivative(Matrix(Q)) : missing
    numerical_rank = !exact && n <= exact_rank_max ? rank(Matrix{Float64}(Q)) : missing
    tested_rank = coalesce(exact_rank,numerical_rank)
    rank_ok = ismissing(tested_rank) || tested_rank == n-1
    reflection_ok = !reflection || (scaled_residual(DR,opnorm(Dhat,Inf)) <= tol &&
                                     scaled_residual(HR,opnorm(Hhat,Inf)) <= tol)
    order_ok = exact_order ? all(==(minimum_order),taylor.orders) :
                             all(>=(minimum_order),taylor.orders)
    expected_B = spdiagm(0 => [i==1 ? -1 : i==n ? 1 : 0 for i in 1:n])
    endpoint_ok = scaled_residual(op.B-expected_B,1) <= tol
    certified = endpoint_ok && sbp_relative <= tol && positive && order_ok && reflection_ok && rank_ok
    return (; certified, exact, positive, symmetric, endpoint_ok, rank_ok, exact_rank, numerical_rank,
            reflection_ok, order_ok, sbp_absolute, sbp_relative, tol, residual=E,
            D_reflection=DR,H_reflection=HR,taylor,Dhat,Hhat,Q)
end

function norm_diagnostics(op::OperatorData; degree=16)
    n = length(op.x)
    Hhat = op.H ./ op.h
    values = isdiag(Hhat) ? sort(Float64.(diag(Hhat))) : eigvals(Symmetric(Matrix{Float64}(Hhat)))
    # Normalize coordinates to [0, 1] to avoid artificial scaling with origin.
    y = (op.x .- first(op.x)) ./ (last(op.x)-first(op.x))
    X = hcat((y.^k for k in 0:degree)...)
    scale = last(op.x)-first(op.x)
    bilinear = X' * (op.H * X) - [scale/(i+j+1) for i in 0:degree,j in 0:degree]
    m = op.depth
    blocks = m > 0 && isdiag(Hhat[m+1:n-m,m+1:n-m]) &&
             all(iszero,Hhat[1:m,m+1:n]) && all(iszero,Hhat[1:n-m,n-m+1:n])
    bounds = blocks ? extrema(vcat(eigvals(Symmetric(Matrix{Float64}(Hhat[1:m,1:m]))),
        eigvals(Symmetric(Matrix{Float64}(Hhat[n-m+1:n,n-m+1:n]))),
        Float64.(diag(Hhat)[m+1:n-m]))) : missing
    return (; eigenvalues=values, min=first(values),max=last(values),
             condition=last(values)/first(values), weights=diag(Hhat),
             quadrature=vec(bilinear[1,:]), bilinear, block_bounds=bounds)
end

function stencil_diagnostics(Dhat, Hhat, Q)
    n = size(Dhat,1)
    radius, counts = zeros(Int,n), zeros(Int,n)
    maxcoeff, l1, l2, distance, signed = (zeros(n) for _ in 1:5)
    I,J,V = findnz(Dhat)
    for (i,j,value) in zip(I,J,V)
        a, offset = Float64(abs(value)), j-i
        radius[i] = max(radius[i],abs(offset))
        counts[i] += 1
        maxcoeff[i] = max(maxcoeff[i],a)
        l1[i] += a
        l2[i] = hypot(l2[i],a)
        distance[i] += a*abs(offset)
        signed[i] += a*offset
    end
    denominator = ifelse.(iszero.(l1),1,l1)
    return (; radius, counts, maxcoeff,l1,l2, displacement=distance./denominator,
            signed_displacement=signed./denominator, D=bandwidth(Dhat),
            H=bandwidth(Hhat), Q=bandwidth(Q), qmax=maximum(abs,Q), nnz=nnz(Dhat),
            cancellation=l1./max.(maxcoeff,eps(Float64)))
end

function row_groups(op)
    n,m = length(op.x),op.depth
    bulk = m+1:n-m
    return (; all=1:n,left=1:m,right=n-m+1:n,bulk,
        transition=unique(vcat(first(bulk,min(2,length(bulk))),last(bulk,min(2,length(bulk))))),
        phase_a=filter(isodd,bulk),phase_b=filter(iseven,bulk))
end

function differentiation_errors(op::OperatorData; precision=256, wavenumbers=(1,4,8))
    return setprecision(precision) do
        T = eltype(op.D) <: Rational ? BigFloat : Float64
        D,H,x = SparseMatrixCSC{T,Int}(op.D),SparseMatrixCSC{T,Int}(op.H),T.(op.x)
        samples = [("exp", exp.(x), exp.(x)),
                   ("chirp",cos.(x+x.^2),-(1 .+ 2x).*sin.(x+x.^2))]
        gaussian = exp.(-20 .* (x .- T(0.5)).^2)
        push!(samples,("gaussian",gaussian,-40 .* (x .- T(0.5)).*gaussian))
        for k in wavenumbers
            push!(samples,("sin_$k",sinpi.(2k .* x),2k*T(π).*cospi.(2k .* x)))
        end
        return map(samples) do (name,u,ux)
            e = D*u-ux
            groups = map(ids -> isempty(ids) ? missing : maximum(abs,e[ids]),row_groups(op))
            (; name,error=Float64.(e), weighted=Float64(sqrt(real(dot(e,H*e)))),groups...)
        end
    end
end

function diagnose_operator(op::OperatorData; minimum_order=1, degree=28,
                           exact_order=false, reflection=true, exact_rank_max=129,
                           spectrum=true, sat_penalty=1.0, numerical_rank_rtol=nothing)
    algebra = algebraic_diagnostics(op;minimum_order,degree,exact_order,reflection,exact_rank_max)
    algebra.certified || return (; algebra)
    norms = norm_diagnostics(op;degree=degree÷2)
    stencil = stencil_diagnostics(algebra.Dhat,algebra.Hhat,algebra.Q)
    kernel = spectrum ? nullspace_diagnostics(op;rtol=numerical_rank_rtol) : nothing
    spectral = spectrum ? spectral_diagnostics(op;sat_penalty) : nothing
    return (; algebra,norms,stencil,kernel,spectral)
end

diagnose_operator(op; depth=0,period=1,kwargs...) =
    diagnose_operator(OperatorData(op;depth,period);kwargs...)
