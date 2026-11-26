# Render one operator at a time so Cairo and FreeType caches are released
# when the worker exits. Completed operator jobs are durable on disk.
function main(args)
    length(args)==1 || (length(args)==2 && args[2]=="--force") ||
        error("usage: generate_figures.jl output_dir [--force]")
    out=abspath(first(args))
    force=length(args)==2
    isfile(joinpath(out,"campaign.jls")) || error("complete the diagnostic stage before rendering")
    dir=joinpath(out,"figures")
    mkpath(dir)
    operators=filter(id->isfile(joinpath(out,"operators",id,"figure_data.jls")),
                     readdir(joinpath(out,"operators")))
    julia=Base.julia_cmd()
    script=joinpath(@__DIR__,"render_figures.jl")
    for id in vcat(operators,["_comparison"])
        marker=joinpath(dir,".$id.complete")
        !force && isfile(marker) && continue
        force && rm(marker;force=true)
        println("Rendering $id")
        flush(stdout)
        run(`$julia --project=$(@__DIR__) --startup-file=no --threads=1 --heap-size-hint=2G $script $out $id`)
        touch(marker)
    end
    println("Figures saved to $dir")
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main(ARGS)
