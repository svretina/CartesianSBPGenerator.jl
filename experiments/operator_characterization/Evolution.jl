function endpoints(op)
    n = length(op.x)
    H = Matrix{Float64}(op.H)
    E = zeros(n,2)
    E[1,1] = E[n,2] = 1
    lifts = Symmetric(H)\E
    return (; H,E,lifts,x=Float64.(op.x),n)
end

function advection_model(op; penalty=1.0, ppw=nothing)
    (; H,E,lifts,x,n) = endpoints(op)
    D = Matrix{Float64}(op.D)
    L = -D - penalty*lifts[:,1]*E[:,1]'
    if isnothing(ppw)
        modes = 2π .* [1,3]
        # The auxiliary oscillators encode the exact time-dependent inflow.
        M = zeros(n+4,n+4)
        M[1:n,1:n] .= L
        M[1:n,n+1] .= penalty*lifts[:,1]
        M[1:n,n+4] .= 0.25penalty*lifts[:,1]
        for (j,w) in enumerate(modes)
            a,b = n+2j-1,n+2j
            M[a,b],M[b,a] = -w,w
        end
        reference = t -> vcat(sinpi.(2 .* (x.-t)) .+ 0.25cospi.(6 .* (x.-t)),
                            -sin(modes[1]*t),cos(modes[1]*t),
                            -sin(modes[2]*t),cos(modes[2]*t))
        initial = reference(0.0)
        kind = "inflow"
    else
        k = 2π/(Float64(op.h)*ppw)
        center = first(x) + 0.32*(last(x)-first(x))
        halfwidth = 0.28*(last(x)-first(x))
        shape(y) = abs(y)<1 ? exp(-4y^2-y^4/(1-y^2)) : 0.0
        packet(y) = shape((y-center)/halfwidth)*cos(k*(y-center))
        reference = t -> packet.(x.-t)
        initial = reference(0.0)
        M = L
        kind = "packet"
    end
    inflow(u) = isnothing(ppw) ? u[n+1]+0.25u[n+4] : 0.0
    energy(u) = real(dot(view(u,1:n),H*view(u,1:n)))/2
    expected_power(u,t) = (abs2(u[1])-abs2(u[n]))/2 -
                           penalty*real(conj(u[1])*(u[1]-inflow(u)))
    mass(u) = sum(H*view(u,1:n))
    mass_power(u,t) = u[1]-u[n]-penalty*(u[1]-inflow(u))
    power(u,rhs) = real(dot(view(u,1:n),H*view(rhs,1:n)))
    return (;kind,M,initial,reference,n,H,energy,expected_power,power,mass,mass_power,
             report_field=1:n, parameters=(;penalty,ppw))
end

function wave_model(op; mode=1, absorbing=false)
    (; H,E,lifts,x,n) = endpoints(op)
    D=Matrix{Float64}(op.D)
    right,left = n+1:2n,2n+1:3n
    M = zeros(3n,3n)
    M[1:n,right] .= I(n)/2
    M[1:n,left] .= I(n)/2
    M[right,right] .= -D - lifts[:,1]*E[:,1]'
    M[left,left] .= D - lifts[:,2]*E[:,2]'
    if !absorbing
        M[right,left] .= -lifts[:,1]*E[:,1]'
        M[left,right] .= -lifts[:,2]*E[:,2]'
    end
    function reference(t)
        if absorbing
            y = (x.-0.32.-t)./0.25
            phi = [abs(q)<1 ? exp(1-inv(1-q*q)) : 0.0 for q in y]
            psi = [abs(q)<1 ? -8q*phi[i]/(1-q*q)^2 : 0.0 for (i,q) in enumerate(y)]
            return vcat(phi,-2psi,zeros(n))
        end
        phi = sinpi.(mode.*x).*cospi(mode*t)
        pi_field = -mode*π .* sinpi.(mode.*x).*sinpi(mode*t)
        psi = mode*π .* cospi.(mode.*x).*cospi(mode*t)
        return vcat(phi,pi_field-psi,pi_field+psi)
    end
    energy(u)=sum(dot(view(u,r),H*view(u,r)) for r in (right,left))/4
    function expected_power(u,t)
        a,b,c,d = u[first(right)],u[last(right)],u[first(left)],u[last(left)]
        return absorbing ? -(abs2(a)+abs2(b)+abs2(c)+abs2(d))/4 :
                           -(abs2(a+c)+abs2(b+d))/4
    end
    constraint(u) = maximum(abs,(u[left]-u[right])/2-D*view(u,1:n))
    profile = sinpi.(mode.*x)
    function mode_coefficient(u)
        pi_field = (u[left]+u[right])/2
        return complex(dot(profile,H*view(u,1:n)),dot(profile,H*pi_field)/(mode*π))
    end
    power(u,rhs) = sum(dot(view(u,r),H*view(rhs,r)) for r in (right,left))/2
    return (;kind=absorbing ? "wave_absorbing" : "standing_wave", M,
            initial=reference(0.0),reference,n,H,energy,expected_power,power,constraint,mode_coefficient,
            report_field=1:n, parameters=(;mode, penalty=1.0))
end

function variable_model(op; compatible, degenerate=false)
    (;H,x,n) = endpoints(op)
    a = degenerate ? 1 .- x : 1 .+ 0.5sinpi.(2x)
    ap = degenerate ? fill(-1.0,n) : π.*cospi.(2x)
    shape = sinpi.(2x)
    force = -shape + 2π.*a.*cospi.(2x) + 0.5ap.*shape
    split = variable_split(op,a;compatible)
    M = [split.L (force+a[1]*shape[1]*split.lift); zeros(1,n) -1.0]
    reference(t) = exp(-t).*vcat(shape,1.0)
    energy(u) = dot(view(u,1:n),H*view(u,1:n))/2
    expected_power(u,t) = (a[1]*abs2(u[1])-a[end]*abs2(u[n]))/2 -
                          a[1]*u[1]*(u[1]-u[end]*shape[1]) +
                          dot(view(u,1:n),H*force)*u[end]
    power(u,rhs) = dot(view(u,1:n),H*view(rhs,1:n))
    return (;kind=compatible ? "variable_compatible" : "variable_standard",M,
            initial=reference(0.0),reference,n,H,energy,expected_power,power,report_field=1:n,
            parameters=(;compatible,degenerate,split.commutator,split.residual,split.abscissa))
end

# The periodic grid has no duplicated endpoint and wraps complete stencil cells.
function periodic_model(op,taps; mode=1, pure=false)
    n = length(op.x)-1
    p = length(taps)
    iszero(n%p) || throw(ArgumentError("periodic grid must contain complete stencil periods"))
    h=Float64(op.h)
    x=Float64.(op.x[1:n].-first(op.x))
    omega = 2π*mode/(n*h)
    function wrap(key)
        indices = [(i,mod1(i+offset,n),Float64(value)) for i in 1:n
                    for (offset,value) in getproperty(taps[mod1(i,p)],key)]
        return sparse(first.(indices),getindex.(indices,2),last.(indices),n,n)
    end
    W,Q=wrap(:norm),wrap(:derivative)
    H=h.*W
    F=cholesky(Symmetric(W))
    M=-(F\Matrix(Q))/h
    basis = cis.(omega.*x)
    if pure && p==2
        Wtheta,Qtheta=Diagnostics.symbol(taps,h*omega)
        U=cholesky(Hermitian(Wtheta)).U
        eig=eigen(Hermitian(-im*(U'\Qtheta/U)))
        idx=argmax(abs2.(eig.vectors'*(U*ones(p))))
        v=U\eig.vectors[:,idx]
        basis .*= [v[mod1(i,p)] for i in 1:n]
    end
    reference(t)=imag.(basis.*cis(-omega*t))
    initial=reference(0.0)
    energy(u)=dot(u,H*u)/2
    expected_power(u,t)=0.0
    mass(u)=sum(H*u)
    mass_power(u,t)=0.0
    power(u,rhs) = dot(u,H*rhs)
    return (;kind=pure ? "periodic_bloch" : "periodic_sampled",M,initial,reference,n,H,energy,
             expected_power,power,mass,mass_power,report_field=1:n,parameters=(;mode,pure,period=p))
end

function rk4_step!(u,A,step,k1,k2,k3,k4,tmp)
    mul!(k1,A,u)
    @. tmp=u+(step/2)*k1
    mul!(k2,A,tmp)
    @. tmp=u+(step/2)*k2
    mul!(k3,A,tmp)
    @. tmp=u+step*k3
    mul!(k4,A,tmp)
    @. u+=(step/6)*(k1+2k2+2k3+k4)
    return nothing
end

struct UnstableEvolution <: Exception
    abscissa::Float64
    final_time::Float64
end
Base.showerror(io::IO,e::UnstableEvolution) =
    print(io,"spectral abscissa ",e.abscissa," predicts amplification beyond exp(20) by T=",e.final_time)

function evolve(model,times;max_steps=200000,safety=0.7,cfl=nothing,
                temporal_order=2,accuracy_factor=0.1,temporal_fraction=0.05,refinements=3)
    (;M,initial,n,reference,report_field)=model
    a = sparse(M)
    values = isnothing(cfl) ? eigvals(M) : nothing
    step_limit = isnothing(cfl) ? rk4_limit(values) : cfl
    h = get(model.parameters,:h,1/max(n-1,1))
    # Shrink RK4 steps for convergence experiments above order four.
    accuracy_step = accuracy_factor*h^max(1,temporal_order/4)
    timestep = min(safety*step_limit,accuracy_step)
    intervals = diff(times)
    demand = timestep > 0 ? sum(ceil.(intervals./timestep)) : Inf
    explicit = 0 < demand <= max_steps
    u=copy(initial)
    scratch=ntuple(_ -> similar(u),5)
    explicit && rk4_step!(copy(initial),a,zero(timestep),scratch...)
    states = zeros(length(u),length(times))
    states[:,1] .= u
    runtimes = zeros(length(times))
    # Distinct output intervals can require distinct exponential propagators.
    propagators = Dict{Float64,Matrix{Float64}}()
    for (i,dt) in enumerate(intervals)
        if explicit
            count = ceil(Int,dt/timestep)
            runtimes[i+1] = runtimes[i]+@elapsed for _ in 1:count
                rk4_step!(u,a,dt/count,scratch...)
            end
        else
            if !isnothing(values) && maximum(real,values)*last(times)>20
                throw(UnstableEvolution(maximum(real,values),last(times)))
            end
            propagator=get!(propagators,dt) do
                exp(dt.*M)
            end
            u=propagator*u
        end
        all(isfinite,u) || error("evolution produced nonfinite state")
        states[:,i+1] .= u
    end
    # An independently computed semidiscrete solution quantifies temporal error.
    exact_states = reduce(hcat, reference.(times))
    errors = states-exact_states
    reference_final = exp(last(times).*M)*initial
    temporal_error = explicit ? maximum(abs,states[report_field,end]-reference_final[report_field]) : missing
    half=exp((last(times)/2).*M)
    exponential_defect=maximum(abs,half*(half*initial)-reference_final)
    spatial_error = maximum(abs,reference_final[report_field]-exact_states[report_field,end])
    temporal_target=max(temporal_fraction*spatial_error,1e-12*max(1,maximum(abs,reference_final[report_field])))
    time_resolved=!explicit || temporal_error<=temporal_target
    if !time_resolved && refinements>0 && 2demand+length(intervals)<=max_steps
        reduced_factor=timestep/(2h^max(1,temporal_order/4))
        return evolve(model,times;max_steps,safety,cfl=step_limit,temporal_order,
            accuracy_factor=reduced_factor,temporal_fraction,refinements=refinements-1)
    end
    energy = map(model.energy,eachcol(states))
    power = [model.expected_power(v,t) for (v,t) in zip(eachcol(states),times)]
    measured_power = [model.power(v,M*v) for v in eachcol(states)]
    integral = vcat(0.0,cumsum(intervals.*(power[1:end-1]+power[2:end])/2))
    energy_balance = energy.-first(energy).-integral
    mass = hasproperty(model,:mass) ? map(model.mass,eachcol(states)) : fill(missing,length(times))
    mass_residual = hasproperty(model,:mass) ?
        [model.mass(M*v)-model.mass_power(v,t) for (v,t) in zip(eachcol(states),times)] :
        fill(missing,length(times))
    constraint = hasproperty(model,:constraint) ? map(model.constraint,eachcol(states)) :
                 fill(missing,length(times))
    fields=states[report_field,:]
    field_errors=errors[report_field,:]
    period = get(model.parameters,:period,0)
    x = collect((0:n-1)./(period>0 ? n : n-1))
    test_mode = cis.(2π.*x)
    projection(u) = dot(test_mode,model.H*view(u,report_field))
    coefficient = hasproperty(model,:mode_coefficient) ? model.mode_coefficient : projection
    phase = [angle(coefficient(view(states,:,i))*conj(coefficient(view(exact_states,:,i))))
             for i in axes(fields,2)]
    for i in 2:length(phase)
        phase[i] += 2π*round((phase[i-1]-phase[i])/(2π))
    end
    weighted = [sqrt(max(dot(e,model.H*e),0)) for e in eachcol(field_errors)]
    maxerror = vec(maximum(abs,field_errors;dims=1))
    parity = (-1.0).^(0:n-1)
    gridscale = [abs(dot(parity,model.H*e))/sqrt(dot(parity,model.H*parity))
                 for e in eachcol(field_errors)]
    pairs = n÷2
    phase_a = field_errors[1:2:2pairs,:]
    phase_b = field_errors[2:2:2pairs,:]
    pair_mean = vec(maximum(abs,(phase_a+phase_b)/2;dims=1))
    pair_difference = vec(maximum(abs,(phase_a-phase_b)/2;dims=1))
    return (;times,x,period,states=fields,reference=exact_states[report_field,:],errors=field_errors, weighted,maxerror,
            max_all_state=vec(maximum(abs,errors;dims=1)),
            step_limit,timestep,requested_steps=demand, explicit, runtimes,temporal_error,spatial_error,
            time_resolved,temporal_target,
            energy,power,power_residual=measured_power-power,energy_balance,mass,mass_residual,
            constraint,phase,gridscale,pair_mean,pair_difference,exponential_defect)
end

function branch_content(values,taps;window=true)
    p = length(taps)
    n = length(values)-mod(length(values),p)
    cells=n÷p
    chi = window ? sinpi.((0:n-1)./max(1,n-1)).^4 : ones(n)
    spectra = fft(reshape(values[1:n].*chi,p,cells),2)./cells
    rows = NamedTuple[]
    for j in 0:cells-1
        k = j <= cells÷2 ? j : j-cells
        theta = 2π*k/n
        H,Q = Diagnostics.symbol(taps,theta)
        U = cholesky(Hermitian(H)).U
        F = eigen(Hermitian(-im*(U'\Q/U)))
        Hp,Qp=Diagnostics.symbol(taps,theta,true)
        v = spectra[:,j+1].*cis.(-theta .* (0:p-1))
        energy = abs2.(F.vectors'*(U*v))/p
        primary = argmax(abs2.(F.vectors'*(U*ones(p))))
        for branch in 1:p
            z=U\F.vectors[:,branch]
            group=real(dot(z,(-im*Qp-F.values[branch]*Hp)*z))
            push!(rows,(; theta,branch=branch==primary ? "physical" : "folded",
                          energy=energy[branch],group_velocity=group))
        end
    end
    return rows
end
