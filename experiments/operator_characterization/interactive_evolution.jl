# Optional desktop viewer. Run with -i in an environment containing GLMakie.
using GLMakie, Serialization

function evolution_index(root)
    entries=Dict{String,Tuple{String,Int}}()
    for (dir,_,files) in walkdir(joinpath(root,"operators"))
        "evolution.jls" in files || continue
        label=replace(relpath(dir,joinpath(root,"operators")), '/' => ' ')
        cache=joinpath(dir,"evolution.jls")
        for (i,record) in enumerate(deserialize(cache))
            entries["$label  •  $(record.tag)"]=(cache,i)
        end
    end
    isempty(entries) && error("no saved evolution data in $root")
    return entries
end

function evolution_viewer(root)
    entries=evolution_index(root)
    labels=sort!(collect(keys(entries)))
    fig=Figure(size=(950,430))
    menu=Menu(fig[1,1:2];options=labels,default=first(labels),width=400)
    slider=Slider(fig[3,1:2];range=0:0.001:1,startvalue=0)
    evolution=lift(menu.selection) do label
        cache,index=entries[label]
        return deserialize(cache)[index].result
    end
    selected=lift(evolution,slider.value) do e,fraction
        idx=argmin(abs.(e.times .- fraction*last(e.times)))
        return (;e,idx)
    end
    Label(fig[4,1:2],lift(z->"Time $(round(z.e.times[z.idx];sigdigits=6))",selected))
    ax=Axis(fig[2,1];xlabel="Position x",ylabel="Field",title="Numerical and analytical solutions")
    error_axis=Axis(fig[2,2];xlabel="Position x",ylabel="Numerical − analytical",title="Pointwise error")
    function coordinates(key)
        return lift(selected) do z
            Point2f.(z.e.x,getproperty(z.e,key)[:,z.idx])
        end
    end
    lines!(ax,coordinates(:reference);color=:black,linestyle=:dash,label="Analytical")
    lines!(ax,coordinates(:states);color=:steelblue,label="Numerical")
    lines!(error_axis,coordinates(:errors);color=:darkorange)
    axislegend(ax)
    on(selected) do _
        autolimits!(ax)
        autolimits!(error_axis)
    end
    return fig
end

length(ARGS)==1 || error("usage: interactive_evolution.jl output_dir")
display(evolution_viewer(only(ARGS)))
