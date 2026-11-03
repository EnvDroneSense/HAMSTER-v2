module HAMSTER

using Dates
using Interpolations
using MPI
using NCDatasets
using Statistics
using TOML

# Land-atmosphere moisture and heat tracking. Modules are included below in the order they are used.

"Package version, from Project.toml."
version() = pkgversion(@__MODULE__)

include("utilities.jl")
include("io.jl")
include("preprocessing.jl")
include("tracking.jl")

export midpoints, cells, only_theta, ongrid, ongrid1d
export read_config, read_mask, partoutput_files, read_window, read_obs, read_observations
export delta, diagnose, correct_deltas, grid_definition, gridindex
export attribute, parcel_blocks, shift_window, track, write_attribution
export write_diagnosis, write_corrected, preprocess_window, preprocess

end
