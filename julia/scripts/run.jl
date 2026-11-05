# Start script: julia --project=julia julia/scripts/run.jl <preprocess|track> <namelist.toml> [blocks]
# Tracking over MPI ranks: mpiexec -n <ranks> julia --project=julia julia/scripts/run.jl track <namelist.toml> [blocks]
using HAMSTER, MPI

function main(args)
    if !(length(args) in (2, 3)) || !(args[1] in ("preprocess", "track"))
        println("usage: julia --project=julia julia/scripts/run.jl <preprocess|track> <namelist.toml> [blocks]")
        return 1 # wrong input
    end
    config = read_config(args[2])
    println("HAMSTER ", HAMSTER.version(), ": ", args[1], " with ", args[2])
    start = time()
    rank = 0
    if args[1] == "preprocess"
        preprocess(config)
    else
        MPI.Init()
        rank = MPI.Comm_rank(MPI.COMM_WORLD)
        track(config; blocks=length(args) == 3 ? parse(Int, args[3]) : 1)
    end
    @info "Finished $(args[1]) (rank $rank): total time $(round(time() - start; digits=1)) s, peak memory $(round(Sys.maxrss() / 2^20; digits=1)) MiB"
    return 0
end

exit(main(ARGS))
