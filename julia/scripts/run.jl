# Start script: julia --project=julia julia/scripts/run.jl <preprocess|track> <namelist.input>
using HAMSTER

function main(args)
    if length(args) != 2 || !(args[1] in ("preprocess", "track"))
        println("usage: julia --project=julia julia/scripts/run.jl <preprocess|track> <namelist.input>")
        return 1
    end
    println("HAMSTER ", HAMSTER.version(), ": ", args[1], " with ", args[2], " (not implemented yet)")
    return 0
end

exit(main(ARGS))
