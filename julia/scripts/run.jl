# Start script: julia --project=julia julia/scripts/run.jl <preprocess|track> <namelist.toml> [blocks]
# Tracking over MPI ranks: mpiexec -n <ranks> julia --project=julia julia/scripts/run.jl track <namelist.toml> [blocks]
using HAMSTER, MPI

function main(args)
    if !(length(args) in (2, 3)) || !(args[1] in ("preprocess", "track"))
        println("usage: julia --project=julia julia/scripts/run.jl <preprocess|track> <namelist.toml> [blocks]")
        return 1
    end
    config = read_config(args[2])
    println("HAMSTER ", HAMSTER.version(), ": ", args[1], " with ", args[2])
    if args[1] == "preprocess"
        preprocess(config)
    else
        MPI.Init()
        track(config; blocks=length(args) == 3 ? parse(Int, args[3]) : 1)
    end
    return 0
end

exit(main(ARGS))
