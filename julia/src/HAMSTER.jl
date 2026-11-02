module HAMSTER

# Julia port of the Python code in ../python. Modules are added here as they are ported.

"Package version, from Project.toml."
version() = pkgversion(@__MODULE__)

include("utilities.jl")

export midpoints, only_theta, ongrid, ongrid1d, hmin

end
