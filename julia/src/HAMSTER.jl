module HAMSTER

using Dates
using Interpolations
using NCDatasets
using Statistics

# Julia port of the Python code in ../python. Modules are added here as they are ported.

"Package version, from Project.toml."
version() = pkgversion(@__MODULE__)

include("utilities.jl")
include("io.jl")
include("preprocessing.jl")
include("tracking.jl")

export midpoints, only_theta, ongrid, ongrid1d, hmin
export read_mask, partoutput_files, read_window, read_obs, read_observations
export delta, diagnose, correct_deltas, grid_definition, gridindex
export attribute, track_window, track, write_attribution
export write_diagnosis, write_corrected, preprocess_window, preprocess

end
