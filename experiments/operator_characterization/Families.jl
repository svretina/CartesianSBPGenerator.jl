const DIRECT = Dict(2 => (5,2,17), 4 => (9,4,25), 6 => (15,6,37))
const COMPACT = Dict(4 => (6,3,2,17), 6 => (10,4,3,25))
const RF_ORDERS = Dict(2=>4,4=>6,6=>8)
const R = Rational{BigInt}
const CFG_FAMILIES = ("standard","restricted_full","alternating","analytical","compact")
const E = CartesianSBPGenerator.Experimental

exact_decimal(x) = rationalize(BigInt,x)
function amplitude(text)
    parts = parse.(BigInt,split(text,'/'))
    length(parts)==1 && return R(only(parts))
    length(parts)==2 || throw(ArgumentError("invalid amplitude multiplier: $text"))
    return parts[1]//parts[2]
end

function construct_cases(cfg)
    selection = cfg["operators"]["families"]
    orders = cfg["operators"]["orders"]
    cases, failures = NamedTuple[],NamedTuple[]
    function build(factory, family, id, target, interior, multiplier=missing)
        start = time_ns()
        try
            generator,depth,period,stages = factory()
            push!(cases,(; id,family,target,interior,multiplier,depth,period,generator,
                          stages,construction_s=(time_ns()-start)*1e-9))
        catch ex
            ex isa InterruptException && rethrow()
            push!(failures,(; id,n=0,status="construction_failed",reason=sprint(showerror,ex)))
        end
    end
    function conventional(interior,boundary,norm)
        t = @timed design_sbp(interior,boundary;norm,skew_reconstruction=:analytical)
        family = t.value
        n=2*family.closure_width+1
        template=OperatorData(instantiate(family;N=n);depth=family.closure_depth)
        stages=[(;stage="closure_search_and_reconstruction",arithmetic="exact",
                 seconds=t.time,allocated_bytes=t.bytes,status="success")]
        return (n,domain) -> transfer_local(template,n,domain), family.closure_depth, 1, stages
    end
    function seeded(p,period,depth,radius,band,seed_n)
        assembly = @timed E.assemble_stencil_experiment(p;N=seed_n,period,closure_depth=depth,
                                                        stencil_radius=radius,norm_bandwidth=band)
        system = assembly.value
        selection = @timed E.search_stencil(system)
        result = selection.value
        result.status===:success || error("seed status $(result.status)")
        stages=[(;stage,arithmetic="exact",seconds=t.time,allocated_bytes=t.bytes,status="success")
                for (stage,t) in (("seed_assembly",assembly),("seed_solve_and_selection",selection))]
        if band==0
            template=OperatorData(result.operator;depth,period)
            return (n,domain) -> transfer_local(template,n,domain), depth,period,stages
        end
        return (n,domain) -> E.refine_stencil(system,result.operator;N=n,domain), depth,period,stages
    end
    if "standard" in selection
        for (i,b) in unique(vcat([(2p,p) for p in orders],[(2,1),(4,2),(6,3)]))
            build(() -> conventional(i,b,:diagonal),"STD","STD-$i-$b",b,i)
        end
    end
    for p in orders
        if "restricted_full" in selection
            i = RF_ORDERS[p]
            build(() -> conventional(i,p,:restricted_full),"RF","RF-$i-$p",p,i)
        end
        if "alternating" in selection
            depth,radius,n = DIRECT[p]
            build(() -> seeded(p,2,depth,radius,0,n),"ALT","ALT-$p",p,p)
        end
        if "analytical" in selection
            for (i,m) in enumerate(cfg["operators"]["amplitude_multipliers"])
                r = amplitude(m)
                build("ANA","ANA-$p-$i",p,p,r) do
                    t = @timed design_sbp(2p,p;norm=:diagonal,skew_reconstruction=:analytical)
                    baseline = t.value
                    ep = r/(12*big(2)^(2p))
                    depth=max(baseline.closure_depth,p+2)
                    n=max(2*baseline.closure_width+1,2depth+3)
                    seed=E.analytical_perturbation(baseline;N=n,epsilon=ep)
                    template=OperatorData(seed;depth,period=2)
                    stages=[(;stage="baseline_construction",arithmetic="exact",seconds=t.time,
                             allocated_bytes=t.bytes,status="success")]
                    return (n,domain) -> transfer_local(template,n,domain),depth,2,stages
                end
            end
        end
        if "compact" in selection
            if haskey(COMPACT,p)
                depth,radius,band,n=COMPACT[p]
                build(() -> seeded(p,1,depth,radius,band,n),"CMP","CMP-$p",p,p)
            else
                push!(failures,(; id="CMP-$p",n=0,status="unavailable",reason="no existing compact seed for this order"))
            end
        end
        push!(failures,(; id="ALT-B-$p",n=0,status="excluded",reason="banded direct construction was dropped from scope"))
    end
    return cases,failures
end

function profile_conventional(case,op)
    types=(R,Float64)
    norm=case.family=="STD" ? :diagonal : :restricted_full
    spec=SBPSpec(interior_order=case.interior,boundary_order=case.target,
                 norm,closure_depth=case.depth)
    rows=NamedTuple[]
    for (arithmetic,T) in zip((:exact,:numerical),types)
        assemble_constraints(spec;arithmetic,method=:reduced)
        build=@timed assemble_constraints(spec;arithmetic,method=:reduced)
        system=build.value
        solve_constraints(system)
        solve=@timed solve_constraints(system)
        block=Matrix{T}(op.H[1:case.depth,1:case.depth]./op.h)
        positivity=@timed CartesianSBPGenerator._positive_definite(block)
        skew=@timed reconstruct_skew(system,block;method=:analytical)
        for (stage,t,status) in (("fixed_depth_assembly",build,"success"),
            ("fixed_depth_linear_solve",solve,string(solve.value.status)),
            ("norm_positivity",positivity,string(positivity.value)),
            ("skew_reconstruction",skew,string(skew.value.status)))
            push!(rows,(;stage,arithmetic=string(arithmetic),seconds=t.time,
                         allocated_bytes=t.bytes,status))
        end
    end
    return rows
end

# Transfer the local repeating Ĥ/Q template, avoiding dense exact rank
# computations on large grids. The diagnostic gate independently checks
# the assembled matrices and records the exact-rank resolution cutoff.
function transfer_local(seed::OperatorData, n, domain)
    n0 = length(seed.x)
    m,p = seed.depth,seed.period
    n >= n0 && mod(n-n0,p)==0 || throw(ArgumentError("refinement must preserve the seed phases"))
    h = (last(domain)-first(domain))/(n-1)
    x = collect(first(domain):h:last(domain))
    H0 = seed.H ./ seed.h
    Q0 = dropzeros!(seed.H*seed.D)
    function extend(A)
        I,J,V = Int[],Int[],eltype(A)[]
        for i in 1:n
            s = i <= m ? i : i > n-m ? n0+i-n : m+mod1(i-m,p)
            ids,values = findnz(A[s,:])
            for (j,value) in zip(ids,values)
                push!(I,i); push!(J,i+j-s); push!(V,value)
            end
        end
        return sparse(I,J,V,n,n)
    end
    H,Q=extend(H0),extend(Q0)
    if isdiag(H)
        D = spdiagm(0 => inv.(diag(H))./h)*Q
    else
        bulk=m+1:n-m
        left,right=1:m,n-m+1:n
        isdiag(H[bulk,bulk]) && iszero(nnz(H[left,m+1:n])) &&
            iszero(nnz(H[1:n-m,right])) ||
            throw(ArgumentError("compact refinements use the exact norm solve"))
        D = spdiagm(0 => inv.(diag(H))./h)*Q
        for rows in (left,right)
            columns=sort!(unique(findnz(Q[rows,:])[2]))
            D[rows,columns] = sparse(Matrix(H[rows,rows])\Matrix(Q[rows,columns])/h)
        end
        dropzeros!(D)
    end
    return OperatorData(D,h.*H,spdiagm(0=>[i==1 ? -1 : i==n ? 1 : 0 for i in 1:n]),x;
                        depth=m,period=p)
end
