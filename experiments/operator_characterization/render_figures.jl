using Serialization, TOML

include("Figures.jl")

length(ARGS)==2 || error("usage: render_figures.jl output_dir operator_id")
out,id=ARGS
config=TOML.parsefile(joinpath(out,"config.toml"))
Figures.generate(out,config,id)
println("Renderer peak RSS for $id: $(round(Sys.maxrss()/2^20;digits=1)) MiB")
