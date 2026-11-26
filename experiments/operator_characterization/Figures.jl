"""Static figures regenerated from persisted benchmark data."""
module Figures

using CairoMakie, LinearAlgebra, Serialization, Statistics, FFTW

const COLORS = Dict("STD"=>:black,"RF"=>"#0072B2","ALT"=>"#D55E00",
                    "ANA"=>"#009E73","CMP"=>"#CC79A7")
const PHASE = ("#0072B2","#D55E00","#009E73","#CC79A7")
const INK = "#202530"
const GRID = "#E7EAEE"
const FLOOR = 1e-30

function publication_theme(font)
    variant(style) = CairoMakie.Makie.assetpath("fonts","TeXGyreHerosMakie-$style.otf")
    fonts = font=="TeX Gyre Heros" ?
        (;regular=variant("Regular"),bold=variant("Bold"),italic=variant("Italic")) :
        (;regular=font,bold=font*" Bold",italic=font*" Italic")
    return Theme(;fonts,
        fontsize=10,figure_padding=(12,12,10,12),backgroundcolor=:white,
        textcolor=INK,Axis=(;titlesize=10,titlefont=:bold,xlabelsize=10,ylabelsize=10,
            xticklabelsize=8,yticklabelsize=8,xgridcolor=GRID,ygridcolor=GRID,
            spinewidth=0.65,xtickwidth=0.65,ytickwidth=0.65,xticksize=3,yticksize=3,
            topspinevisible=false,rightspinevisible=false),
        Colorbar=(;labelsize=9,ticklabelsize=8,width=9),
        Legend=(;labelsize=8,framevisible=false,patchsize=(18,8)),
        Lines=(;linewidth=1.35))
end

function export_figure(fig,dir,name,opts)
    opts["save_pdf"] && save(joinpath(dir,"$name.pdf"),fig;pt_per_unit=1)
    opts["save_png"] && save(joinpath(dir,"$name.png"),fig;px_per_unit=opts["png_dpi"]/72)
    GC.gc()
end

positive(x) = max(abs(Float64(x)),FLOOR)
signed_log(x,scale) = sign(x)*log10(1+abs(x)/scale)
newfig(title;size=(510,410)) = begin
    fig=Figure(;size)
    Label(fig[0,1:2],title;font=:bold,fontsize=11,halign=:left,tellwidth=false)
    fig
end
function panel(fig,r,c,title,x,y;kwargs...)
    return Axis(fig[r,c];title,xlabel=x,ylabel=y,kwargs...)
end
function map_panel(fig,r,c,title,z;symbol="",matrix=true,mode=:signed,scale=0.01)
    grid=fig[r,c]=GridLayout()
    ax=Axis(grid[1,1];title, xlabel=matrix ? "Column" : "Degree",
            ylabel=matrix ? "Row" : "Node",yreversed=true)
    A=permutedims(Float64.(z))
    if mode===:pattern
        heatmap!(ax,Float64.(.!iszero.(A));colormap=[:white,INK],colorrange=(0,1),rasterize=5)
    else
        logvalues=mode===:absolute ? log10.(positive.(A)) : signed_log.(A,scale)
        if mode===:absolute
            lo,hi=extrema(logvalues)
            colorrange=(min(lo,hi-1),hi)
            cmap=:viridis
        else
            peak=max(maximum(abs,logvalues),0.1)
            colorrange=(-peak,peak)
            cmap=:balance
        end
        image=heatmap!(ax,logvalues;colormap=cmap,colorrange,rasterize=5)
        label=mode===:absolute ? "log₁₀ |$symbol|" : "sgn($symbol) log₁₀(1 + |$symbol|/$scale)"
        Colorbar(grid[1,2],image;label)
    end
    return ax
end

function operator_label(r)
    id = "$(r.id)  ($(r.bulk_orders), $(r.minimum_order))"
    return ismissing(r.amplitude) ? id : "$id, ε=$(replace(string(r.amplitude),"//"=>"/"))"
end

function operator_figures(data,dir,opts)
    r=data.summary
    title="$(operator_label(r))  •  N=$(r.n)"
    x=data.x
    a,b=extrema(x)
    col=COLORS[r.family]
    fig=newfig(title;size=(510,465))
    map_panel(fig,1,1,"Derivative sparsity",data.D;mode=:pattern)
    map_panel(fig,1,2,"Signed derivative coefficients",data.D;symbol="hD")
    map_panel(fig,2,1,"Norm sparsity",data.H;mode=:pattern)
    map_panel(fig,2,2,"Weighted derivative sparsity",data.Q;mode=:pattern)
    ax=panel(fig,3,1,"Dimensionless norm weights","Position x","diag(H)/h")
    lines!(ax,x,Float64.([v.weight for v in data.structure]);color=col)
    ax=panel(fig,3,2,"Selected stencils","Offset from evaluation node","hD coefficient")
    for (j,i) in enumerate(unique([1,min(5,r.n),cld(r.n,2),cld(r.n,2)+1]))
        js=findall(!iszero,data.D[i,:])
        scatterlines!(ax,js.-i,data.D[i,js];color=PHASE[j],markersize=3,label="Row $i")
    end
    axislegend(ax;position=:rt)
    export_figure(fig,dir,"structure",opts)
    if r.family=="CMP"
        fig=newfig(title;size=(510,380))
        map_panel(fig,1,1,"Compact norm",data.H;symbol="H/h")
        map_panel(fig,1,2,"Weighted derivative",data.Q;symbol="Q")
        map_panel(fig,2,1,"Explicit derivative",data.D;symbol="hD")
        ax=panel(fig,2,2,"Norm eigenvalues","Index","λ(H/h)")
        lines!(ax,data.eigenvalues;color=col)
        export_figure(fig,dir,"compact_matrices",opts)
    end

    for (name,A,subject,symbol) in (("sbp_residual",data.sbp,"SBP identity residual","R"),
                                   ("reflection_residual",data.Dreflection,"Derivative reflection residual","Rᴰ"))
        fig=newfig(title;size=(510,220))
        map_panel(fig,1,1,subject,A;symbol,mode=:absolute)
        if name=="reflection_residual"
            map_panel(fig,1,2,"Norm reflection residual",data.Hreflection;symbol="Rᴴ",mode=:absolute)
        else
            ax=panel(fig,1,2,"Norm eigenvalue distribution","Index","λ(H/h)")
            lines!(ax,data.eigenvalues;color=col)
        end
        export_figure(fig,dir,name,opts)
    end

    fig=newfig(title;size=(510,475))
    ax=panel(fig,1,1,"Pointwise accuracy","Position x","Order")
    stairs!(ax,x,[v.order for v in data.structure];color=col)
    ax=panel(fig,1,2,"Row-wise coefficient norms","Position x","Norm";yscale=log10)
    for (j,(key,label)) in enumerate(((:coefficient_max,"Maximum"),(:coefficient_l1,"ℓ₁"),(:coefficient_l2,"ℓ₂")))
        lines!(ax,x,positive.(getproperty.(data.structure,key));color=PHASE[j],label)
    end
    axislegend(ax;position=:rt)
    ax=panel(fig,2,1,"Stencil geometry","Position x","Grid offsets")
    lines!(ax,x,[v.stencil_radius for v in data.structure];color=PHASE[1],label="Radius")
    lines!(ax,x,[v.displacement for v in data.structure];color=PHASE[2],label="Mean distance")
    lines!(ax,x,[v.signed_displacement for v in data.structure];color=PHASE[3],label="Signed displacement")
    axislegend(ax;position=:ct)
    ax=panel(fig,2,2,"Leading Taylor coefficient","Position x","Signed coefficient")
    lead=[ismissing(v.leading) ? NaN : Float64(v.leading) for v in data.structure]
    for (j,phase) in enumerate((1:2:length(x),2:2:length(x)))
        lines!(ax,x[phase],lead[phase];color=PHASE[j],label=j==1 ? "Phase A" : "Phase B")
    end
    axislegend(ax;position=:rt)
    ax=panel(fig,3,1,"Magnitude of leading coefficient","Position x","|Leading coefficient|";yscale=log10)
    lines!(ax,x,positive.(lead);color=col)
    ax=panel(fig,3,2,"Smooth differentiation error","Position x","D f − f′")
    for (j,item) in enumerate(first(data.errors,min(3,length(data.errors))))
        lines!(ax,x,item.error;color=PHASE[j],label=item.name)
    end
    axislegend(ax;position=:rt)
    export_figure(fig,dir,"row_accuracy",opts)

    fig=newfig(title;size=(510,352))
    map_panel(fig,1,1,"Taylor residuals through degree $(size(data.moments,2)-1)",data.moments;
              matrix=false,mode=:absolute,symbol="C")
    ax=panel(fig,1,2,"Scalar quadrature residual","Polynomial degree","Absolute error";yscale=log10)
    scatterlines!(ax,0:length(data.qres)-1,positive.(data.qres);color=col,markersize=4)
    ax=panel(fig,2,1,"Maximum Taylor residual","Polynomial degree","max |C|";yscale=log10)
    scatterlines!(ax,0:size(data.moments,2)-1,
                  positive.(vec(maximum(abs,data.moments;dims=1)));color=col,markersize=3)
    map_panel(fig,2,2,"Bilinear norm residual",data.qbilinear;symbol="E",mode=:absolute)
    export_figure(fig,dir,"polynomial_residuals",opts)
    return
end

function spectrum_figures(data,dir,opts)
    isnothing(data.specs) && return
    r=data.summary
    spec=data.specs[argmin(abs.([s.sat_penalty for s in data.specs].-1))]
    k=data.kernel
    title="$(operator_label(r))  •  N=$(r.n)  •  τ=$(spec.sat_penalty)"
    fig=newfig(title;size=(510,495))
    ax=panel(fig,1,1,"Raw derivative","Re λ(hD)","Im λ(hD)")
    scatter!(ax,real.(spec.raw),imag.(spec.raw);color=PHASE[1],markersize=4)
    ax=panel(fig,1,2,"Inflow SAT operator","Re λ(hL)","Im λ(hL)")
    numerical_range=get(spec,:numerical_range,ComplexF64[])
    if !isempty(numerical_range)
        lines!(ax,real.(numerical_range),imag.(numerical_range);
               color=:gray65,linestyle=:dash,label="Energy numerical range")
    end
    scatter!(ax,real.(spec.sat),imag.(spec.sat);color=PHASE[2],markersize=4)
    isempty(numerical_range) || axislegend(ax;position=:lt)
    ax=panel(fig,2,1,"RK4 stability region","Re z","Im z")
    xx,yy=range(-3.2,0.5;length=151),range(-3.3,3.3;length=151)
    rk(z)=abs(1+z*(1+z*(1/2+z*(1/6+z/24))))
    contour!(ax,xx,yy,[rk(complex(x,y)) for x in xx,y in yy];levels=[1.0],color=:gray40)
    scatter!(ax,real.(spec.sat).*spec.cfl,imag.(spec.sat).*spec.cfl;color=PHASE[2],markersize=4)
    ax=panel(fig,2,2,"SAT penalty sensitivity","Penalty τ","Dimensionless CFL")
    scatterlines!(ax,[s.sat_penalty for s in data.specs],[s.cfl for s in data.specs];color=PHASE[2],label="Spectral")
    if !ismissing(spec.energy_cfl)
        scatterlines!(ax,[s.sat_penalty for s in data.specs],[s.energy_cfl for s in data.specs];
                      color=PHASE[1],label="Energy, sampled")
    end
    axislegend(ax;position=:rt)
    ax=panel(fig,3,1,"Most boundary-localized SAT mode","Position x","|Eigenvector|")
    lines!(ax,data.x,abs.(spec.most_localized);color=PHASE[2])
    ax=panel(fig,3,2,"Homogeneous amplification","Dimensionless time t/h","Induced energy norm")
    hlines!(ax,[1.0];color=:gray70,linestyle=:dash)
    if !isempty(spec.transient)
        lines!(ax,collect(spec.times),spec.transient;color=PHASE[1],label="Semidiscrete")
        scatter!(ax,collect(spec.steps).*spec.safety.*spec.cfl,spec.amplification;
                 color=PHASE[2],markersize=5,label="RK4")
        axislegend(ax;position=:rt)
    end
    export_figure(fig,dir,"spectrum",opts)

    fig=newfig(title;size=(510,350))
    ax=panel(fig,1,1,"Singular-value spectra","Index","Singular value";yscale=log10)
    lines!(ax,k.singular .|> positive;color=PHASE[1],label="Dimensionless")
    lines!(ax,k.energy_singular .|> positive;color=PHASE[2],label="Energy transform")
    axislegend(ax;position=:lb)
    ax=panel(fig,1,2,"Near-null singular vectors","Position x","Amplitude")
    for j in axes(k.near_vectors,2)
        lines!(ax,data.x,k.near_vectors[:,j];color=PHASE[j])
    end
    ax=panel(fig,2,1,"Near-null Fourier power","Discrete Fourier index","Power";yscale=log10)
    powers=abs2.(fft(k.near_vectors,1))
    for j in axes(powers,2)
        lines!(ax,0:r.n-1,positive.(powers[:,j]);color=PHASE[j])
    end
    ax=panel(fig,2,2,"Checkerboard differentiation","Position x","D v")
    lines!(ax,data.x,k.checkerboard_derivative[:,1];color=PHASE[1],label="Full")
    lines!(ax,data.x,k.checkerboard_derivative[:,2];color=PHASE[2],label="Windowed")
    axislegend(ax;position=:rt)
    export_figure(fig,dir,"nullspace",opts)
    if !isnothing(data.ps)
        ps=data.ps
        fig=newfig(title;size=(510,300))
        layout=fig[1,1:2]=GridLayout()
        ax=Axis(layout[1,1];title="Energy-norm pseudospectrum",xlabel="Re z",ylabel="Im z")
        Z=log10.(positive.(ps.values))
        field=heatmap!(ax,ps.x,ps.y,Z;colormap=:viridis,rasterize=5)
        contour!(ax,ps.x,ps.y,Z;levels=Float64.(ps.levels),color=:white,linewidth=0.75)
        scatter!(ax,real.(spec.sat),imag.(spec.sat);color=:white,markersize=3)
        Colorbar(layout[1,2],field;label="log₁₀ σₘᵢₙ(zI − A)")
        export_figure(fig,dir,"pseudospectrum",opts)
    end
end

function dispersion_figure(data,dir,opts)
    isnothing(data.dispersion) && return
    d=data.dispersion.data
    r=data.summary
    title="$(operator_label(r))  •  N=$(r.n)"
    fig=newfig(title;size=(510,485))
    for (i,(key,label,units)) in enumerate(((:frequencies,"Modified wavenumber","κ"),
            (:phase,"Phase velocity","cₚ / c"),(:group,"Group velocity","cᵍ / c"),
            (:branch_energy,"Smooth Fourier-mode projection","Energy fraction")))
        logarithmic = i==4 && size(d.branch_energy,2)>1
        ax=panel(fig,cld(i,2),mod1(i,2),label,"Wavenumber θ","$units";
                 yscale=logarithmic ? log10 : identity)
        values=getproperty(d,key)
        i==1 && lines!(ax,d.thetas,d.thetas;color=:gray65,linestyle=:dash)
        i in (2,3) && hlines!(ax,[1.0];color=:gray65,linestyle=:dash)
        for j in axes(values,2)
            series = logarithmic ? positive.(values[:,j]) : values[:,j]
            lines!(ax,d.thetas,series;color=PHASE[j],label=j==1 ? "Physical" : "Folded")
        end
        i==1 && axislegend(ax;position=:lt)
    end
    ax=panel(fig,3,1,"Low-wavenumber phase error","θ","|κ/θ − 1|";xscale=log10,yscale=log10)
    low=filter(x->x.phase_error>0,data.dispersion.low)
    scatterlines!(ax,getproperty.(low,:theta),getproperty.(low,:phase_error);
                  color=PHASE[1],markersize=3)
    ax=panel(fig,3,2,"Resolution thresholds","Relative error tolerance","Points per wavelength";xscale=log10)
    thresholds=collect(d.thresholds)
    for (j,(ppw,label)) in enumerate(((d.ppw_phase,"Phase"),(d.ppw_group,"Group")))
        indices=findall(isfinite,collect(ppw))
        scatterlines!(ax,thresholds[indices],collect(ppw)[indices];color=PHASE[j],label,markersize=4)
    end
    axislegend(ax;position=:rt)
    export_figure(fig,dir,"dispersion",opts)
end

function evolution_figures(data,records,dir,opts)
    r=data.summary
    for (tag,kind,e) in records
        # Heatmaps have nonuniform time coordinates when requested output times are inserted.
        x=e.x
        title="$(operator_label(r))  •  N=$(r.n)  •  $tag"
        fig=newfig(title;size=(510,486))
        ax=panel(fig,1,1,"Final numerical and analytical profiles","Position x","Field")
        lines!(ax,x,e.reference[:,end];color=:black,linestyle=:dash,label="Exact")
        lines!(ax,x,e.states[:,end];color=PHASE[1],label=e.explicit ? "RK4" : "Semidiscrete")
        axislegend(ax;position=:rt)
        for (row,col,Z,name) in ((1,2,e.errors,"Space–time error"),(2,1,e.states,"Space–time solution"))
            grid=fig[row,col]=GridLayout()
            ax=Axis(grid[1,1];title=name,xlabel="Position x",ylabel="Time")
            vmax=max(maximum(abs,Z),eps(Float64))
            map=heatmap!(ax,x,e.times,Z;colormap=:balance,colorrange=(-vmax,vmax),rasterize=5)
            Colorbar(grid[1,2],map)
        end
        ax=panel(fig,2,2,"Energy balance","Time","Normalized quantity")
        normalization=max(e.energy[1],eps(Float64))
        lines!(ax,e.times,e.energy./normalization;color=PHASE[1],label="Energy")
        lines!(ax,e.times,e.energy_balance./normalization;color=PHASE[2],label="Balance defect")
        axislegend(ax;position=:rt)
        ax=panel(fig,3,1,"Instantaneous power residual","Time","Absolute residual";yscale=log10)
        lines!(ax,e.times,positive.(e.power_residual);color=PHASE[1])
        ax=panel(fig,3,2,"Phase and grid-scale error","Time","Magnitude";yscale=log10)
        lines!(ax,e.times,positive.(e.phase);color=PHASE[1],label="Phase (rad)")
        lines!(ax,e.times,positive.(e.gridscale);color=PHASE[2],label="Grid-scale projection")
        if !all(ismissing,e.constraint)
            lines!(ax,e.times,positive.(e.constraint);color=PHASE[3],label="Wave constraint")
        end
        axislegend(ax;position=:rt)
        export_figure(fig,dir,"evolution_$tag",opts)
    end
end

function comparison_series!(ax,rows,x,y;label=true,positive_y=false,
                            markersize=4,line=:solid,by=:id)
    ids=unique(getproperty.(rows,by))
    for id in ids
        r=sort(filter(z->getproperty(z,by)==id && !ismissing(getproperty(z,y)),rows);
               by=z->getproperty(z,x))
        isempty(r) && continue
        xs=getproperty.(r,x)
        ys=getproperty.(r,y)
        finite=findall(i->isfinite(xs[i]) && isfinite(ys[i]),eachindex(xs))
        isempty(finite) && continue
        args=label ? (;label=string(id)) : (;)
        series=findfirst(==(first(r).family),("STD","RF","ALT","ANA","CMP"))
        styles=(:solid,:dash,:dot,:dashdot,:dashdotdot)
        markers=(:circle,:rect,:utriangle,:diamond,:cross)
        # Always use the same family color. Distinguish amplitudes by stroke and symbol.
        suffix=tryparse(Int,last(split(string(id),'-')))
        style=first(r).family=="ANA" && !isnothing(suffix) ? mod1(suffix+1,5) : series
        scatterlines!(ax,xs[finite],positive_y ? positive.(ys[finite]) : ys[finite];
            color=COLORS[first(r).family],linestyle=line===:solid ? styles[style] : line,
            marker=markers[style],markersize,args...)
    end
end

function comparison_figures(results,dir,opts)
    for target in sort!(unique(getproperty.(results.summaries,:target)))
        s=filter(r->r.target==target&&isodd(r.n),results.summaries)
        d=filter(r->r.target==target&&isodd(r.n)&&r.function_name=="exp",results.differentiated)
        e=filter(r->r.target==target&&isodd(r.n),results.evolution)
        isempty(s) && continue
        title="Equal minimum pointwise accuracy  •  Order $target"
        where=joinpath(dir,"order_$target")
        mkpath(where)
        fig=newfig(title;size=(510,550))
        for (i,(key,label,unit)) in enumerate(((:norm_condition,"Norm condition number","κ₂(H/h)"),
                (:sat_radius,"SAT spectral radius","ρ(hL)"),(:cfl,"RK4 spectral step limit","CFL"),
                (:physical_sigma_nonzero,"Smallest nonzero singular value","σₙ₋₁(D)"),
                (:abte,"Average boundary Taylor coefficient","ABTE"),
                (:mbte,"Maximum boundary Taylor coefficient","MBTE")))
            ax=panel(fig,cld(i,2),mod1(i,2),label,"Nodes N",unit;xscale=log10,yscale=log10)
            comparison_series!(ax,s,:n,key;positive_y=true,label=i==1)
            i==1 && axislegend(ax;position=:rt,nbanks=2)
        end
        export_figure(fig,where,"operator_comparison",opts)

        fig=newfig(title;size=(510,340))
        for (i,(key,label,unit)) in enumerate(((:norm_min,"Smallest norm eigenvalue","λₘᵢₙ(H/h)"),
                  (:norm_max,"Largest norm eigenvalue","λₘₐₓ(H/h)"),
                  (:max_coeff,"Explicit derivative coefficients","Maximum |hD|"),
                  (:max_Q,"Weighted derivative coefficients","Maximum |Q|")))
            ax=panel(fig,cld(i,2),mod1(i,2),label,"Nodes N",unit;xscale=log10,yscale=log10)
            comparison_series!(ax,s,:n,key;positive_y=true,label=i==1)
            i==1 && axislegend(ax;position=:rt,nbanks=2)
        end
        export_figure(fig,where,"norm_and_coefficients",opts)

        resolved=filter(r->!ismissing(r.sat_radius) && !ismissing(r.bulk_radius),s)
        if !isempty(resolved)
            rows=map(r->merge(r,(;radius_ratio=r.sat_radius/r.bulk_radius)),resolved)
            fig=newfig(title;size=(510,245))
            ax=panel(fig,1,1:2,"Finite-domain spectral radius relative to the bulk",
                     "Nodes N","ρ(hL) / max |κ|";xscale=log10,yscale=log10)
            hlines!(ax,[1.0];color=:gray65,linestyle=:dash)
            comparison_series!(ax,rows,:n,:radius_ratio;positive_y=true)
            axislegend(ax;position=:rt,nbanks=2)
            export_figure(fig,where,"finite_domain_vs_bulk",opts)
        end

        fig=newfig(title;size=(510,495))
        for (i,region) in enumerate((:left,:bulk,:phase_a,:phase_b,:weighted,:transition))
            ax=panel(fig,cld(i,2),mod1(i,2),"$(uppercasefirst(replace(string(region),'_'=>' '))) error",
                     "Grid spacing h","Differentiation error";xscale=log10,yscale=log10)
            rows=filter(r->r.region==region,d)
            comparison_series!(ax,rows,:h,:error;positive_y=true,label=i==1)
            i==1 && !isempty(rows) && axislegend(ax;position=:lt,nbanks=2)
            if !isempty(rows)
                sample=last(sort(rows;by=r->r.h))
                xx=extrema([r.h for r in rows])
                slopes=region in (:left,:transition,:weighted) ? (target,) : (target,2target)
                for (j,p) in enumerate(slopes)
                    lines!(ax,collect(xx),sample.error.*(collect(xx)./sample.h).^p;
                           color=:gray55,linestyle=j==1 ? :dash : :dot,linewidth=0.8)
                end
            end
        end
        export_figure(fig,where,"differentiation_convergence",opts)

        fig=newfig(title;size=(510,320))
        ax=panel(fig,1,1,"Fixed-resolution application cost","Nodes N","Best measured application (s)";xscale=log10,yscale=log10)
        comparison_series!(ax,s,:n,:apply_s;positive_y=true)
        axislegend(ax;position=:lt,nbanks=2)
        driven=map(r -> merge(r,(;apply_s=only(filter(q->q.id==r.id&&q.n==r.n,s)).apply_s)),filter(z->z.region==:all,d))
        ax=panel(fig,1,2,"Differentiation work–precision","Application time (s)","Maximum error";xscale=log10,yscale=log10)
        comparison_series!(ax,driven,:apply_s,:error;positive_y=true,label=false)
        export_figure(fig,where,"differentiation_work_precision",opts)
        isempty(e) && continue
        variants=("inflow","periodic_sampled","wave_absorbing","standing_wave","variable_standard","variable_compatible")
        for (kind,metric,label) in ((:convergence,:h,"Grid spacing h"),(:work_precision,:runtime_s,"Measured RK4 runtime (s)"))
            fig=newfig(title;size=(510,505))
            for (i,problem) in enumerate(variants)
                ax=panel(fig,cld(i,2),mod1(i,2),replace(uppercasefirst(problem),'_'=>' '),label,
                         "Maximum error";xscale=log10,yscale=log10)
                rows=filter(r->r.problem==problem &&
                    (kind==:convergence || (r.rk4 && get(r,:time_resolved,true))),e)
                isempty(rows) && continue
                primary=first(sort!(unique(getproperty.(rows,:tag))))
                rows=filter(r->r.tag==primary,rows)
                tf=maximum(getproperty.(rows,:time))
                rows=filter(r->r.time==tf,rows)
                ax.title="$(ax.title[])  •  T=$tf"
                comparison_series!(ax,rows,metric,:error;positive_y=true,label=i==1)
                i==1 && axislegend(ax;position=:lt,nbanks=2)
            end
            export_figure(fig,where,"pde_$kind",opts)
        end
        packet=filter(r->r.problem=="packet",e)
        if !isempty(packet)
            finest=maximum(getproperty.(packet,:n))
            packet=filter(r->r.n==finest && r.time==maximum(getproperty.(packet,:time)),packet)
            rows=map(r->merge(r,(;ppw=parse(Int,last(split(r.tag,'-'))))),packet)
            fig=newfig("$title  •  N=$finest";size=(510,340))
            ax=panel(fig,1,1:2,"Energy remaining after the packet exits",
                     "Carrier points per wavelength","E(T) / E(0)";yscale=log10)
            comparison_series!(ax,rows,:ppw,:energy;positive_y=true)
            axislegend(ax;position=:rt,nbanks=2)
            for (i,(key,name)) in enumerate(((:reflected_physical,"Physical branch"),(:reflected_folded,"Folded branch")))
                ax=panel(fig,2,i,name,"Carrier points per wavelength","Incoming / incident energy";
                         yscale=log10)
                comparison_series!(ax,rows,:ppw,key;positive_y=true,label=false)
            end
            export_figure(fig,where,"boundary_reflection",opts)
        end
    end
end

function generate(out,cfg,id::AbstractString)
    options=cfg["output"]
    CairoMakie.activate!()
    dir=joinpath(out,"figures")
    mkpath(dir)
    with_theme(publication_theme(options["font"])) do
        if id=="_comparison"
            results=deserialize(joinpath(out,"campaign.jls"))
            comparison_figures(results,dir,options)
        else
            folder=joinpath(out,"operators",id)
            data=deserialize(joinpath(folder,"figure_data.jls"))
            dest=joinpath(dir,data.summary.id)
            mkpath(dest)
            operator_figures(data,dest,options)
            spectrum_figures(data,dest,options)
            dispersion_figure(data,dest,options)
            evolution=joinpath(folder,string(data.summary.n),"evolution.jls")
            if isfile(evolution)
                records=deserialize(evolution)
                selected=[r for r in records if !startswith(r.tag,"packet-") || r.tag=="packet-8"]
                evolution_figures(data,selected,dest,options)
            end
        end
    end
    return dir
end

end
