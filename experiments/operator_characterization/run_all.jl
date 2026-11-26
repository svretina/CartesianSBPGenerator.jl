# Keep the diagnostics and Cairo renderer in separate address spaces.
function main(args)
    options=filter(a->!startswith(a,"--"),args)
    flags=filter(a->startswith(a,"--"),args)
    all(in(("--no-pde","--no-plots","--resume")),flags) && length(options)<=2 ||
        error("usage: run_all.jl [config.toml] [output_dir] [--no-pde] [--no-plots] [--resume]")
    config=isempty(options) ? joinpath(@__DIR__,"config.toml") : options[1]
    out=length(options)==2 ? options[2] : joinpath(@__DIR__,"results")
    julia=Base.julia_cmd()
    phase_flags=filter(!=("--no-plots"),flags)
    run(`$julia --project=$(@__DIR__) --startup-file=no --threads=1 --heap-size-hint=2G $(joinpath(@__DIR__,"characterize.jl")) $config $out $phase_flags`)
    if !("--no-plots" in flags)
        run(`$julia --project=$(@__DIR__) --startup-file=no --threads=1 $(joinpath(@__DIR__,"generate_figures.jl")) $out`)
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main(ARGS)
