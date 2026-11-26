using CartesianSBPGenerator
using CartesianSBPGenerator.Diagnostics
using LinearAlgebra, SparseArrays, Statistics, SHA, Serialization, TOML, Printf
using FFTW, BenchmarkTools

include("Families.jl")
include("Evolution.jl")
include("Campaign.jl")

function main(args)
    options=filter(a->!startswith(a,"--"),args)
    flags=filter(a->startswith(a,"--"),args)
    all(in(("--no-pde","--resume")),flags) && length(options)<=2 ||
        error("usage: characterize.jl [config.toml] [output_dir] [--no-pde] [--resume]")
    config=read_config(isempty(options) ? DEFAULT_CONFIG : options[1])
    out=length(options)==2 ? options[2] : joinpath(@__DIR__,"results")
    results=campaign(config,out;pde=!("--no-pde" in flags),resume="--resume" in flags)
    println("Completed $(length(results.summaries)) operator/grid pairs. See $out/status.csv.")
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main(ARGS)
