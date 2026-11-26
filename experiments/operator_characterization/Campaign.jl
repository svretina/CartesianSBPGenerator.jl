const DEFAULT_CONFIG = joinpath(@__DIR__,"config.toml")

function read_config(path)
    config = TOML.parsefile(DEFAULT_CONFIG)
    if abspath(path) != abspath(DEFAULT_CONFIG)
        for (section,data) in TOML.parsefile(path)
            merge!(get!(config,section,Dict{String,Any}()),data)
        end
    end
    cfg = config["grid"]
    length(cfg["domain"])==2 && cfg["domain"] == [0,1] ||
        throw(ArgumentError("the manufactured PDE campaign uses [0, 1]"))
    cfg["sizes"] == sort(unique(cfg["sizes"])) && all(>(3),cfg["sizes"]) ||
        throw(ArgumentError("grid sizes must be increasing and distinct"))
    0 < config["time_integration"]["cfl_safety_factor"] < 1 ||
        throw(ArgumentError("the CFL safety factor must lie in (0, 1)"))
    isfinite(config["time_integration"]["temporal_error_fraction"]) &&
        config["time_integration"]["temporal_error_fraction"] > 0 ||
        throw(ArgumentError("temporal error fraction must be positive"))
    operators=config["operators"]
    all(in(CFG_FAMILIES),operators["families"]) ||
        throw(ArgumentError("unknown operator family in configuration"))
    all(in(keys(DIRECT)),operators["orders"]) ||
        throw(ArgumentError("campaign orders must be selected from 2, 4 and 6"))
    issorted(config["pseudospectrum"]["contour_levels"]) ||
        throw(ArgumentError("pseudospectral contour levels must be in ascending order"))
    return config
end

function csv_field(x)
    ismissing(x) && return ""
    s=string(x)
    return any(c -> c in (',','"','\n','\r'),s) ? "\""*replace(s,"\""=>"\"\"")*"\"" : s
end

function write_rows(path, rows)
    isempty(rows) && return
    open(path,"w") do io
        println(io,join(keys(first(rows)),','))
        for row in rows
            println(io,join(csv_field.(values(row)),','))
        end
    end
end

function write_matrix(path,A; rowkey=:row, colkey=:column, column_offset=0, sparse_data=false)
    indices = sparse_data ? zip(findnz(sparse(A))...) :
        ((i,j,A[i,j]) for j in axes(A,2) for i in axes(A,1))
    open(path,"w") do io
        println(io,"$rowkey,$colkey,value")
        for (i,j,v) in indices
            println(io,join((i,j+column_offset,csv_field(v)),','))
        end
    end
end

function coefficient_hash(algebra)
    io=IOBuffer()
    for (name,matrix) in ((:H,algebra.Hhat),(:Q,algebra.Q))
        println(io,name," ",size(matrix,1))
        for (i,j,value) in zip(findnz(matrix)...)
            println(io,i," ",j," ",value)
        end
    end
    return bytes2hex(sha256(take!(io)))
end

function timing(apply!,u,samples)
    out=similar(u)
    apply!(out,u)
    trial=@benchmark $apply!($out,$u) samples=samples evals=1 seconds=30
    stats=median(trial)
    return (;seconds=stats.time*1e-9, bytes=trial.memory, allocations=trial.allocs,
             points_per_second=length(u)/(stats.time*1e-9))
end

function application_performance(op,algebra,samples,cfl)
    D = SparseMatrixCSC{Float64,Int}(op.D)
    Q = SparseMatrixCSC{Float64,Int}(algebra.Q)
    H = SparseMatrixCSC{Float64,Int}(op.H)
    dense = Matrix(D)
    template = sparse(transpose(D))
    u = sinpi.(2Float64.(op.x))
    start=time_ns()
    F = cholesky(Symmetric(H))
    factor_time=(time_ns()-start)*1e-9
    temp=similar(u)
    stencil! = function(out,input)
        fill!(out,0.0)
        rows,values=rowvals(template),nonzeros(template)
        for i in eachindex(out), k in nzrange(template,i)
            out[i]+=values[k]*input[rows[k]]
        end
        return out
    end
    # H⁻¹Q uses the cached factor and records allocations of the actual solve.
    weighted! = function(out,input)
        mul!(temp,Q,input)
        copyto!(out,F\temp)
        return out
    end
    test=similar(u)
    weighted!(test,u)
    isapprox(test,D*u;rtol=1e-7,atol=1e-7) || error("cached norm application differs from D")
    statistics = map((:dense,:sparse,:stencil,:cached_norm),
                      ((out,v)->mul!(out,dense,v),(out,v)->mul!(out,D,v),stencil!,weighted!),
                      (dense,D,template,(F,H,Q))) do method,apply!,storage
        t=timing(apply!,u,samples)
        (;method,t...,storage_bytes=Base.summarysize(storage),
          factor_seconds=method===:cached_norm ? factor_time : 0.0,
          seconds_per_simulated_time=!ismissing(cfl) && cfl>0 ?
              4*t.seconds/(cfl*Float64(op.h)) : missing)
    end
    return collect(statistics)
end

function low_dispersion(taps)
    return setprecision(192) do
        rows=map(BigFloat(2).^(-collect(3:20))) do theta
            H,Q=Diagnostics.symbol(taps,theta)
            B=-im*(H\Q)
            if size(B,1)==1
                frequency=real(B[1,1])
            else
                midpoint=real(tr(B))/2
                spread=sqrt(max(0,midpoint^2-real(det(B))))
                frequency=midpoint+spread
                abs(frequency-theta)>abs(midpoint-spread-theta) && (frequency=midpoint-spread)
            end
            (;theta=Float64(theta),phase_error=Float64(abs(frequency/theta-1)))
        end
        return rows
    end
end

function analytical_reference_error(taps,case,disp)
    p=case.target
    stencil=centered_stencil(;derivative_order=1,accuracy_order=2p)
    epsilon=Float64(case.multiplier/(12*big(2)^(2p)))
    k0(t)=sum(Float64(c)*sin(j*t) for (j,c) in zip(stencil.offsets,coefficients(stencil)))
    return maximum(enumerate(disp.thetas)) do (i,t)
        a,b=k0(t),k0(t+π)
        beta=(-1)^(p÷2+1)*2.0^(p+2)*cos(t)*sin(t)^(p+1)
        midpoint=(a+b)/2
        spread=hypot((a-b)/2,epsilon*beta)
        maximum(abs,sort(disp.frequencies[i,:])-[midpoint-spread,midpoint+spread])
    end
end

function save_structure(out,op,report)
    (;algebra,norms,stencil)=report
    C=algebra.taylor.coefficients
    rows=[(;node=i,x=op.x[i], order=algebra.taylor.orders[i],
           leading=algebra.taylor.leading[i], weight=norms.weights[i],
           stencil_radius=stencil.radius[i], nonzeros=stencil.counts[i],
           coefficient_max=stencil.maxcoeff[i], coefficient_l1=stencil.l1[i], coefficient_l2=stencil.l2[i],
           displacement=stencil.displacement[i],signed_displacement=stencil.signed_displacement[i],
           cancellation=stencil.cancellation[i]) for i in eachindex(op.x)]
    write_rows(joinpath(out,"rows.csv"),rows)
    for (name,A) in ((:D,algebra.Dhat),(:H,algebra.Hhat),(:Q,algebra.Q),(:B,op.B),
        (:sbp_residual,algebra.residual),(:D_reflection,algebra.D_reflection),
        (:H_reflection,algebra.H_reflection))
        write_matrix(joinpath(out,"$name.csv"),A;sparse_data=true)
    end
    write_matrix(joinpath(out,"taylor.csv"),C;colkey=:degree,column_offset=-1)
    write_matrix(joinpath(out,"norm_bilinear.csv"),norms.bilinear)
    write_rows(joinpath(out,"norm.csv"),
        [(;index=i,value=v) for (i,v) in enumerate(norms.eigenvalues)])
    return rows
end

function save_spectrum(out,op,cfg,specs,kernel)
    evrows=[(;penalty=s.sat_penalty,index=i,real=real(z),imag=imag(z),localization=s.localization[i])
            for s in specs for (i,z) in enumerate(s.sat)]
    first_spec=first(specs)
    write_rows(joinpath(out,"sat_eigenvalues.csv"),evrows)
    write_rows(joinpath(out,"raw_eigenvalues.csv"),
        [(;real=real(z),imag=imag(z)) for z in first_spec.raw])
    write_rows(joinpath(out,"sat_summary.csv"),
        [(;penalty=s.sat_penalty,s.radius,s.abscissa,s.strict,s.tol,s.cfl,s.energy_cfl,s.energy_abscissa,s.departure)
         for s in specs])
    write_rows(joinpath(out,"numerical_range.csv"),
        [(;penalty=s.sat_penalty,index=i,real=real(z),imag=imag(z))
         for s in specs for (i,z) in enumerate(s.numerical_range)])
    write_rows(joinpath(out,"transient.csv"),
       [(;penalty=s.sat_penalty,time=t,amplification=y) for s in specs for (t,y) in zip(s.times,s.transient)])
    write_rows(joinpath(out,"rk4_energy.csv"),
       [(;penalty=s.sat_penalty,step=k,amplification=y) for s in specs for (k,y) in zip(s.steps,s.amplification)])
    write_rows(joinpath(out,"singular_values.csv"),
        [(;index=i,dimensionless=s,energy=kernel.energy_singular[i],physical=s/Float64(op.h))
         for (i,s) in enumerate(kernel.singular)])
    write_rows(joinpath(out,"near_null_rank.csv"),
        [(;rtol=t,rank=k) for (t,k) in zip((1e-6,1e-8,1e-10,1e-12,1e-14,1e-16),kernel.ranks)])
    write_matrix(joinpath(out,"near_null_vectors.csv"),kernel.near_vectors)
    write_matrix(joinpath(out,"near_null_fourier.csv"),abs2.(fft(kernel.near_vectors,1)))
    write_matrix(joinpath(out,"checkerboard_derivative.csv"),kernel.checkerboard_derivative)
    main=specs[argmin(abs.([s.sat_penalty for s in specs].-1))]
    write_rows(joinpath(out,"localized_mode.csv"),
        [(;node=i,real=real(v),imag=imag(v)) for (i,v) in enumerate(main.most_localized)])
    options=cfg["pseudospectrum"]
    if options["enabled"] && length(op.x) in options["grid_sizes"]
        radius=max(main.radius,1e-6)
        M=main.energy_matrix
        x=collect(range(-1.1radius,0.15radius;length=options["sampling_points"]))
        y=collect(range(-1.1radius,1.1radius;length=options["sampling_points"]))
        values=pseudospectrum(M,x,y)
        write_rows(joinpath(out,"pseudospectrum.csv"),
            [(;re=a,im=b,sigma_min=values[i,j]) for (j,b) in enumerate(y) for (i,a) in enumerate(x)])
        return (;x,y,values,levels=options["contour_levels"])
    end
    return nothing
end

function save_dispersion(out,taps,case,cfg)
    options=cfg["dispersion"]
    d=dispersion_diagnostics(taps;samples=options["sampling_points"],
                             thresholds=options["relative_error_thresholds"])
    validation=case.family=="ANA" ? analytical_reference_error(taps,case,d) : missing
    ismissing(validation) || validation <= 5e-12 || error("analytical Bloch formula mismatch")
    low=low_dispersion(taps)
    write_rows(joinpath(out,"dispersion.csv"),
        [(;theta=t,branch=j,wavenumber=d.frequencies[i,j],phase=d.phase[i,j],group=d.group[i,j],
           projection=d.branch_energy[i,j],dissipation=0.0,polarization=d.polarization[i,j])
         for (i,t) in enumerate(d.thetas) for j in axes(d.frequencies,2)])
    write_rows(joinpath(out,"dispersion_low.csv"),low)
    write_rows(joinpath(out,"dispersion_ppw.csv"),
        [(;tolerance=t,phase=d.ppw_phase[i],group=d.ppw_group[i])
         for (i,t) in enumerate(options["relative_error_thresholds"])])
    orders=[log2(low[i+1].phase_error/low[i+2].phase_error) for i in 2:5
            if low[i+1].phase_error>0 && low[i+2].phase_error>0]
    return (;data=d,low,validation,phase_order=isempty(orders) ? missing : median(orders))
end

function save_evolution(out,model,result)
    write_rows(joinpath(out,"trace.csv"),
        [(;time=t,energy=result.energy[i],power=result.power[i],power_residual=result.power_residual[i],
           energy_balance=result.energy_balance[i],mass=result.mass[i],mass_residual=result.mass_residual[i],
           constraint=result.constraint[i],phase=result.phase[i],gridscale=result.gridscale[i],
           cell_mean_error=result.pair_mean[i],cell_difference_error=result.pair_difference[i],
           error_inf=result.maxerror[i],error_H=result.weighted[i],error_all_state=result.max_all_state[i],
           runtime_s=result.explicit ? result.runtimes[i] : missing) for (i,t) in enumerate(result.times)])
    write_rows(joinpath(out,"fields.csv"),
        [(;time=t,x=result.x[j],
           numerical=result.states[j,i],reference=result.reference[j,i],error=result.errors[j,i])
         for (i,t) in enumerate(result.times) for j in 1:model.n])
end

function pde_campaign(dir,op,taps,case,cfg,report_failure)
    tc=cfg["time_integration"]
    ad=cfg["advection"]
    wc=cfg["wave"]
    results=NamedTuple[]
    visuals=NamedTuple[]
    function run(make_model,tag,final_times)
        try
            model=make_model()
            target=joinpath(dir,tag)
            mkpath(target)
            open(joinpath(target,"parameters.toml"),"w") do io
                TOML.print(io,Dict(string(k)=>v for (k,v) in pairs(model.parameters) if !isnothing(v)))
            end
            times=sort!(unique(vcat(collect(range(0.0,Float64(maximum(final_times));length=tc["snapshots"])),
                                    Float64.(final_times))))
            result=evolve(model,times;max_steps=tc["max_steps"],safety=tc["cfl_safety_factor"],
                          temporal_order=case.target,temporal_fraction=tc["temporal_error_fraction"])
            save_evolution(target,model,result)
            reflection=(physical=missing,folded=missing)
            if !isnothing(taps) && model.kind=="packet"
                incident=branch_content(result.states[:,1],taps)
                remaining=branch_content(result.states[:,end],taps)
                scale=sum(r.energy for r in incident if r.group_velocity>0;init=0.0)
                if scale>0
                    incoming(branch) = sum(r.energy for r in remaining
                        if r.group_velocity<0 && r.branch==branch;init=0.0)/scale
                    reflection=(physical=incoming("physical"),folded=incoming("folded"))
                end
            end
            for tf in final_times
                i=searchsortedfirst(times,tf)
                ids=Diagnostics.row_groups(op)
                field=result.errors[:,i]
                edge=filter(<=(length(field)),vcat(ids.left,ids.right))
                bulk=filter(<=(length(field)),ids.bulk)
                push!(results,(;id=case.id,case.family,case.target,n=length(op.x),h=Float64(op.h),
                       problem=model.kind,tag,time=tf,error=result.maxerror[i],weighted=result.weighted[i],
                       boundary=isempty(edge) ? missing : maximum(abs,field[edge]),
                       interior=isempty(bulk) ? missing : maximum(abs,field[bulk]),
                       energy=result.energy[i]/result.energy[1],phase=result.phase[i],
                       maximum_power_residual=maximum(abs,result.power_residual),
                       constraint=result.constraint[i],gridscale=result.gridscale[i],
                       reflected_physical=tf==last(times) ? reflection.physical : missing,
                       reflected_folded=tf==last(times) ? reflection.folded : missing,
                       rk4=result.explicit,requested_steps=result.requested_steps, dt=result.timestep,
                       time_resolved=result.time_resolved,
                       temporal_error=tf==last(times) ? result.temporal_error : missing,
                       spatial_error=tf==last(times) ? result.spatial_error : missing,
                       exponential_defect=tf==last(times) ? result.exponential_defect : missing,
                       runtime_s=result.explicit ? result.runtimes[i] : missing))
            end
            if model.kind in ("inflow","standing_wave","wave_absorbing","packet","variable_compatible","variable_standard","periodic_sampled")
                push!(visuals,(;tag,model.kind,result))
            end
            if !isnothing(taps) && model.kind in ("packet","periodic_sampled","periodic_bloch")
                window=model.kind=="packet"
                indices=window ? unique([1,argmin(abs.(times.-0.5)),length(times)]) : collect(eachindex(times))
                rows=[(;time=times[i],r...) for i in indices
                      for r in branch_content(result.states[:,i],taps;window)]
                write_rows(joinpath(target,"branch_energy.csv"),rows)
                reduced=[(;time=times[i],folded=sum(r.energy for r in rows if r.time==times[i] && r.branch=="folded";init=0.0),
                           incoming=sum(r.energy for r in rows if r.time==times[i] && r.group_velocity<0;init=0.0),
                           total=sum(r.energy for r in rows if r.time==times[i];init=0.0)) for i in indices]
                write_rows(joinpath(target,"branch_totals.csv"),reduced)
            end
        catch ex
            report_failure(tag,ex)
        end
    end
    run(()->advection_model(op),"inflow",ad["final_times"])
    T=last(op.x)-first(op.x)
    for ppw in ad["packet_ppw"]
        run(()->advection_model(op;ppw),"packet-$ppw",[0.55T,1.0T,1.25T])
    end
    if !isnothing(taps)
        for pure in (false,true)
            pure && case.period!=2 && continue
            run(()->periodic_model(op,taps;pure),pure ? "periodic-bloch" : "periodic-sampled",T.*ad["crossings"])
        end
    end
    run(()->wave_model(op;absorbing=true),"wave-absorbing",[0.5,1.0])
    for mode in wc["standing_wave_modes"]
        run(()->wave_model(op;mode),"standing-$mode",2 .* wc["number_of_periods"]./mode)
    end
    for compatible in (false,true), degenerate in (false,true)
        tag="variable-$(compatible ? "compatible" : "standard")-$(degenerate ? "degenerate" : "positive")"
        run(()->variable_model(op;compatible,degenerate),tag,[0.25,1.0])
    end
    return results,visuals
end

function save_differentiation(dir,op,errors,case)
    h=Float64(op.h)
    out=NamedTuple[]
    for e in errors
        write_rows(joinpath(dir,"error_$(e.name).csv"),
            [(;node=i,x=Float64(op.x[i]),error=y) for (i,y) in enumerate(e.error)])
        for region in (:all,:left,:right,:transition,:bulk,:phase_a,:phase_b,:weighted)
            value=getproperty(e,region)
            ismissing(value) || push!(out,(;id=case.id,case.family,case.target,n=length(op.x),h,
                                            function_name=e.name,region,error=value))
        end
    end
    return out
end

function add_rates(rows, group)
    previous=Dict{Any,Tuple{Float64,Float64}}()
    return map(sort(rows;by=r->r.n)) do row
        key=group(row)
        current=(row.h,row.error)
        prior=get(previous,key,nothing)
        rate=isnothing(prior) || !(prior[2]>0 && current[2]>0) ? missing :
            log(prior[2]/current[2])/log(prior[1]/current[1])
        previous[key]=current
        merge(row,(;observed_order=rate))
    end
end

function summary_row(op,case,result,secs,spec,kernel,disp,perf)
    n=length(op.x)
    (;algebra,norms,stencil)=result
    bulk=Diagnostics.row_groups(op).bulk
    edge=collect(1:case.depth)
    C=algebra.taylor.coefficients
    errors=Float64.(C[edge,case.target+2])
    weighted_mean, rms = if case.period==2
        j = first(bulk)
        w = norms.weights[j:j+1]
        c = C[j:j+1,case.target+2]
        (Float64(sum(w.*c)/sum(w)),sqrt(Float64(sum(w.*abs2.(c))/sum(w))))
    else
        (missing,missing)
    end
    ppw_index=isnothing(disp) ? nothing : findfirst(==(0.01),disp.data.thresholds)
    return (;id=case.id,case.family,case.target,case.interior,n,h=Float64(op.h),
        parity=iseven(n) ? "even" : "odd",case.depth,case.period,
        amplitude=ismissing(case.multiplier) ? missing : case.multiplier/(12*big(2)^(2*case.target)),
        exact=algebra.exact,hash=coefficient_hash(algebra),minimum_order=minimum(algebra.taylor.orders),
        bulk_orders=join(sort!(unique(algebra.taylor.orders[bulk])),"/"),
        sbp_absolute=Float64(algebra.sbp_absolute),sbp_relative=Float64(algebra.sbp_relative),
        reflection_D=Float64(opnorm(algebra.D_reflection,Inf)),reflection_H=Float64(opnorm(algebra.H_reflection,Inf)),
        certified=algebra.certified, exact_rank=algebra.exact_rank,
        numerical_rank=isnothing(kernel) ? missing : kernel.numerical_rank, positive=algebra.positive,
        norm_min=norms.min,norm_max=norms.max,norm_condition=norms.condition,
        block_bound_min=ismissing(norms.block_bounds) ? missing : first(norms.block_bounds),
        block_bound_max=ismissing(norms.block_bounds) ? missing : last(norms.block_bounds),
        norm_bandwidth=norms.max>0 ? stencil.H.upper : missing, D_lower=stencil.D.lower,D_upper=stencil.D.upper,
        Q_lower=stencil.Q.lower,Q_upper=stencil.Q.upper,nonzeros=stencil.nnz,
        max_coeff=maximum(stencil.maxcoeff),max_Q=Float64(stencil.qmax), min_weight=Float64(minimum(norms.weights)),
        min_order_boundary=case.target,abte=sqrt(sum(abs2,errors)/max(1,length(errors))),
        mbte=maximum(abs,errors;init=0.0),weighted_mean,rms,
        sigma_nonzero=isnothing(kernel) ? missing : kernel.sigma_nonzero,
        physical_sigma_nonzero=isnothing(kernel) ? missing : kernel.physical_sigma_nonzero,
        effective_condition=isnothing(kernel) ? missing : kernel.condition_nonzero,
        checkerboard=isnothing(kernel) ? missing : kernel.checkerboard,
        windowed_checkerboard=isnothing(kernel) ? missing : kernel.windowed_checkerboard,
        sat_penalty=isnothing(spec) ? missing : spec.sat_penalty,
        sat_radius=isnothing(spec) ? missing : spec.radius,
        raw_radius=isnothing(spec) ? missing : spec.raw_radius,
        cfl=isnothing(spec) ? missing : spec.cfl,
        energy_abscissa=isnothing(spec) ? missing : spec.energy_abscissa,
        bulk_radius=isnothing(disp) ? missing : maximum(abs,disp.data.frequencies),
        norm_symbol_min=isnothing(disp) ? missing : disp.data.norm_min,
        bloch_validation=isnothing(disp) ? missing : disp.validation,
        phase_order=isnothing(disp) ? missing : disp.phase_order,
        bulk_skew_residual=isnothing(disp) ? missing : disp.data.skew_residual,
        bulk_sign_crossings=isnothing(disp) ? missing : join(disp.data.crossings,"/"),
        origin_frequencies=isnothing(disp) ? missing : join(disp.data.origin,"/"),
        ppw_1pct_phase=isnothing(ppw_index) ? missing : disp.data.ppw_phase[ppw_index],
        ppw_1pct_group=isnothing(ppw_index) ? missing : disp.data.ppw_group[ppw_index],
        construction_s=case.construction_s,instantiate_s=secs,
        apply_s=minimum(t.seconds for t in perf),memory_bytes=minimum(t.storage_bytes for t in perf))
end

format_metric(x) = ismissing(x) ? "—" : @sprintf("%.4g",x)

function write_campaign_report(out,summaries,evolution,failures)
    open(joinpath(out,"report.md"),"w") do io
        println(io,"# Cartesian SBP characterization\n")
        println(io,"The tables use the largest odd grid with a SAT spectral result in each order group. ",
                "The order pair gives the measured bulk and minimum pointwise orders. ",
                "CFL uses the stability interval of classical RK4 at the reported SAT penalty.\n")
        for order in sort!(unique(getproperty.(summaries,:target)))
            rows=filter(r->r.target==order && isodd(r.n) && !ismissing(r.cfl),summaries)
            isempty(rows) && continue
            n=maximum(getproperty.(rows,:n))
            println(io,"## Order $order on $n nodes\n")
            println(io,"| Operator | Orders | ε | Half-bandwidth | κ₂(H/h) | Max. abs. hD | ABTE | hρ(L) | RK4 CFL | PPW at 1% phase error |")
            println(io,"| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
            for r in filter(q->q.n==n,rows)
                eps_label=ismissing(r.amplitude) ? "—" : replace(string(r.amplitude),"//"=>"/")
                values=(r.id,"($(r.bulk_orders), $(r.minimum_order))",eps_label,
                        max(r.D_lower,r.D_upper),format_metric(r.norm_condition),format_metric(r.max_coeff),
                        format_metric(r.abte),format_metric(r.sat_radius),format_metric(r.cfl),
                        format_metric(r.ppw_1pct_phase))
                println(io,"| ",join(values," | ")," |")
            end
            println(io)
        end
        if !isempty(evolution)
            explicit=count(r->r.rk4,evolution)
            println(io,"## Time-dependent benchmarks\n")
            println(io,"The campaign recorded $(length(evolution)) final-time measurements. ",
                "$explicit used explicit RK4. The others use semidiscrete propagation and are ",
                "excluded from runtime comparisons. The per-case CSV files contain the temporal ",
                "reference error and matrix-exponential consistency check.\n")
        end
        println(io,"## Coverage and status\n")
        println(io,"| Status | Count |\n| --- | ---: |")
        for status in sort!(unique(getproperty.(failures,:status)))
            println(io,"| $status | $(count(r->r.status==status,failures)) |")
        end
        println(io,"\nSee status.csv for per-case reasons. Exact rank above the configured limit ",
            "and sampled dispersion roots require separate verification. ",
            "The packet branch estimates use windowed spectra and include finite-window leakage.")
    end
end

function checkpoint(path,value)
    mkpath(dirname(path))
    mktemp(dirname(path)) do temporary,io
        serialize(io,value)
        close(io)
        mv(temporary,path;force=true)
    end
end

function campaign(cfg,out; pde=true,resume=false)
    populated=isdir(out) && !isempty(readdir(out))
    config_path=joinpath(out,"config.toml")
    metadata_path=joinpath(out,"metadata.toml")
    aggregate=joinpath(out,"campaign.jls")
    if populated
        resume || error("choose a new output directory or pass --resume")
        isfile(config_path) && TOML.parsefile(config_path)==cfg ||
            error("resume requires the original effective configuration")
        meta=TOML.parsefile(metadata_path)
        meta["pde"]==pde || error("resume requires the original PDE selection")
        meta["julia"]==string(VERSION) || error("resume requires the original Julia version")
        isfile(aggregate) && return deserialize(aggregate)
    end
    mkpath(out)
    BLAS.set_num_threads(1)
    if !populated
        open(io -> TOML.print(io,cfg),config_path,"w")
        open(metadata_path,"w") do io
        project=Base.active_project()
        manifest=joinpath(dirname(project),"Manifest.toml")
        root=abspath(joinpath(@__DIR__,"../.."))
        git(args...)=readchomp(Cmd(vcat(["git","-C",root],collect(args))))
        commit=try git("rev-parse","HEAD") catch; "unknown" end
        dirty=try !isempty(git("status","--porcelain","--untracked-files=no")) catch; true end
        cpus=Sys.cpu_info()
        TOML.print(io,Dict("julia"=>string(VERSION),"cpu"=>Sys.CPU_NAME,"machine"=>Sys.MACHINE,
            "cpu_model"=>isempty(cpus) ? "unknown" : first(cpus).model,
            "julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),"git_commit"=>commit,
            "working_tree_modified"=>dirty,"time_integrator"=>"RK4","domain"=>cfg["grid"]["domain"],"pde"=>pde,
            "manifest_sha256"=>isfile(manifest) ? bytes2hex(sha256(read(manifest))) : "missing"))
        end
    end
    sizes=cfg["grid"]["sizes"]
    sizes=cfg["grid"]["test_even_sizes"] ? sort!(unique(vcat(sizes,sizes.+1))) : sizes
    domain=Tuple(exact_decimal.(cfg["grid"]["domain"]))
    differentiated, evolution, summaries, costs, construction = (NamedTuple[] for _ in 1:5)
    println("Constructing configured operator families")
    flush(stdout)
    cases,failures=construct_cases(cfg)
    for case in cases
        append!(construction,[(;id=case.id,n=0,r...) for r in case.stages])
        operator_dir=joinpath(out,"operators",case.id)
        template_path=joinpath(operator_dir,"template.jls")
        seed=resume && isfile(template_path) ? deserialize(template_path) : nothing
        bulk_taps=nothing
        dispersion=nothing
        cached=resume && isfile(joinpath(operator_dir,"figure_data.jls"))
        profiled=false
        for n in sizes
            issue(status,reason)=push!(failures,(;id=case.id,n,status,reason))
            function capture(tag,ex)
                ex isa Union{InterruptException,OutOfMemoryError} && rethrow(ex)
                reason=sprint(showerror,ex)
                issue(ex isa UnstableEvolution ? "numerical_instability" : "$(tag)_failed","$tag: $reason")
                println(stderr,"$(case.id) N=$n $tag: $reason")
            end
            if case.period==2 && iseven(n)
                issue("unsupported_parity","period-two reflection requires odd finite-domain grids")
                continue
            end
            if case.family=="CMP" && n>cfg["grid"]["compact_max_nodes"]
                issue("resource_limit","dense exact compact derivative above configured resolution cap")
                continue
            end
            dir=joinpath(operator_dir,string(n))
            checked=joinpath(dir,"checkpoint.jls")
            streams=(differentiated,evolution,summaries,costs,construction,failures)
            if resume && isfile(checked)
                for (destination,source) in zip(streams,deserialize(checked))
                    append!(destination,source)
                end
                profiled=true
                println("$(case.id)  N=$n (checkpoint)")
                continue
            end
            lengths=map(length,streams)
            mkpath(dir)
            println("$(case.id)  N=$n")
            flush(stdout)
            try
                start=time_ns()
                op=if !isnothing(seed)
                    transfer_local(seed,n,domain)
                else
                    instance=case.generator(n,domain)
                    instance isa OperatorData ? instance :
                        OperatorData(instance;depth=case.depth,period=case.period)
                end
                if isnothing(seed) && case.period==2
                    seed=op
                    checkpoint(template_path,seed)
                end
                secs=(time_ns()-start)*1e-9
                push!(construction,(;id=case.id,n,stage="instantiation",arithmetic="exact",
                    seconds=secs,allocated_bytes=missing,status="success"))
                polynomial_degree=cfg["accuracy"]["maximum_taylor_degree"]
                report=diagnose_operator(op;minimum_order=case.target,degree=polynomial_degree,
                    exact_order=case.family in ("ALT","ANA","CMP"),
                    exact_rank_max=cfg["accuracy"]["exact_rank_max"],spectrum=false)
                a=report.algebra
                if !a.certified
                    issue("invalid","SBP=$(a.sbp_relative), positive=$(a.positive), order=$(a.order_ok), reflection=$(a.reflection_ok), rank=$(a.exact_rank)")
                    write_matrix(joinpath(dir,"sbp_residual.csv"),a.residual;sparse_data=true)
                    continue
                end
                if ismissing(a.exact_rank)
                    issue("exact_rank_not_tested","grid exceeds configured exact-rank limit")
                end
                structure=save_structure(dir,op,report)
                if !profiled && case.family in ("STD","RF")
                    try
                        timings=profile_conventional(case,op)
                        append!(construction,[(;id=case.id,n,r...) for r in timings])
                    catch ex
                        capture("construction_profile",ex)
                    end
                    profiled=true
                end
                errors=differentiation_errors(op;precision=cfg["accuracy"]["bigfloat_precision"],
                    wavenumbers=cfg["accuracy"]["wavenumbers"])
                append!(differentiated,save_differentiation(dir,op,errors,case))
                qrows=[(;degree=k-1,residual=v) for (k,v) in enumerate(report.norms.quadrature)]
                write_rows(joinpath(dir,"quadrature.csv"),qrows)
                specs,kernel,ps = nothing,nothing,nothing
                if n in cfg["spectrum"]["grid_sizes"] || (iseven(n) && n-1 in cfg["spectrum"]["grid_sizes"])
                    kernel=nullspace_diagnostics(op)
                    if kernel.numerical_rank!=n-1
                        issue("numerical_rank_deficiency","SVD numerical rank $(kernel.numerical_rank) at default relative threshold")
                    end
                    advanced=n in cfg["spectrum"]["advanced_sizes"]
                    specs=[spectral_diagnostics(op;sat_penalty=tau,advanced,
                            safety=cfg["time_integration"]["cfl_safety_factor"]) for tau in cfg["spectrum"]["sat_penalties"]]
                    ps=save_spectrum(dir,op,cfg,specs,kernel)
                end
                main=isnothing(specs) ? nothing : specs[argmin(abs.([s.sat_penalty for s in specs].-1))]
                if isnothing(bulk_taps)
                    try
                        taps=bulk_stencil(op)
                        dispersion=save_dispersion(dir,taps,case,cfg)
                        bulk_taps=taps
                    catch ex
                        capture("dispersion",ex)
                    end
                end
                visual=NamedTuple[]
                if pde && n in cfg["time_integration"]["grid_sizes"] &&
                    (case.family!="ANA" || !cfg["time_integration"]["primary_amplitude_only"] || case.multiplier==1)
                    rows,visual=pde_campaign(dir,op,bulk_taps,case,cfg,capture)
                    append!(evolution,rows)
                    for row in rows
                        row.rk4 || issue("rk4_unavailable","$(row.tag) at T=$(row.time): semidiscrete fallback; $(row.requested_steps) requested RK4 steps")
                        row.time_resolved || issue("temporal_error_unresolved","$(row.tag) at T=$(row.time): temporal error target was not reached")
                    end
                    checkpoint(joinpath(dir,"evolution.jls"),visual)
                end
                perf=application_performance(op,a,cfg["performance"]["samples"],isnothing(main) ? missing : main.cfl)
                write_rows(joinpath(dir,"application.csv"),perf)
                append!(costs,[(;id=case.id,n,r...) for r in perf])
                row=summary_row(op,case,report,secs,main,kernel,dispersion,perf)
                push!(summaries,row)
                if !cached || (!isnothing(specs) && n==first(cfg["spectrum"]["grid_sizes"]))
                    artifact=(;summary=row,x=Float64.(op.x),structure,
                              D=Matrix{Float64}(a.Dhat),H=Matrix{Float64}(a.Hhat),Q=Matrix{Float64}(a.Q),
                              sbp=Matrix{Float64}(a.residual),Dreflection=Matrix{Float64}(a.D_reflection),
                              Hreflection=Matrix{Float64}(a.H_reflection),
                              moments=Float64.(a.taylor.coefficients),
                              qres=Float64.(report.norms.quadrature),qbilinear=Float64.(report.norms.bilinear),
                              eigenvalues=report.norms.eigenvalues,
                              errors,specs,kernel,ps,dispersion)
                    written=@timed checkpoint(joinpath(operator_dir,"figure_data.jls"),artifact)
                    push!(construction,(;id=case.id,n,stage="serialization",arithmetic="exact_and_float",
                         seconds=written.time,allocated_bytes=written.bytes,status="success"))
                    cached=true
                end
            catch ex
                capture("diagnostic",ex)
            finally
                rows=map((stream,start)->stream[start+1:end],streams,lengths)
                if !isempty(rows[3]) || any(r->r.status=="invalid",rows[6])
                    checkpoint(checked,rows)
                end
            end
            GC.gc()
        end
    end
    differentiated=add_rates(differentiated,r->(r.id,r.function_name,r.region,iseven(r.n)))
    evolution=add_rates(evolution,r->(r.id,r.tag,r.time,iseven(r.n)))
    write_rows(joinpath(out,"summary.csv"),summaries)
    write_rows(joinpath(out,"differentiation.csv"),differentiated)
    write_rows(joinpath(out,"evolution.csv"),evolution)
    write_rows(joinpath(out,"performance.csv"),costs)
    write_rows(joinpath(out,"construction.csv"),construction)
    write_rows(joinpath(out,"status.csv"),failures)
    for p in sort!(unique([r.target for r in summaries]))
        rows=filter(r->r.target==p,summaries)
        write_rows(joinpath(out,"order_$p.csv"),rows)
    end
    comparisons=NamedTuple[]
    for a in summaries
        a.family=="STD" || continue
        for family in ("RF","ALT","ANA","CMP")
            choices=filter(r->r.n==a.n && r.family==family,summaries)
            isempty(choices) && continue
            for metric in (:nonzeros,:apply_s)
                b=choices[argmin(abs.(log.(getproperty.(choices,metric)./getproperty(a,metric))))]
                push!(comparisons,(;reference=a.id,candidate=b.id,n=a.n,metric,
                      ratio=getproperty(b,metric)/getproperty(a,metric),
                      reference_order=a.minimum_order,candidate_order=b.minimum_order))
            end
        end
    end
    write_rows(joinpath(out,"comparable_cost.csv"),comparisons)
    write_campaign_report(out,summaries,evolution,failures)
    isempty(summaries) && error("no operator/grid pairs completed. See status.csv")
    checkpoint(aggregate,(;summaries,differentiated,evolution,costs,failures))
    println("Diagnostic peak RSS: $(round(Sys.maxrss()/2^20;digits=1)) MiB")
    return (;summaries,differentiated,evolution,costs,failures)
end
