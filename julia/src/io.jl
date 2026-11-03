"Format of the dates in the namelist."
const DATEFORMAT = dateformat"yyyy-mm-dd_HH"

"Time stamp of a file name, `yyyymmddHHMMSS`."
stamp(time) = Dates.format(time, "yyyymmddHHMMSS")

"Keys of the namelist, by section, and the type each value must have."
const NAMELIST = (
    "dates" => ("start_date" => String, "end_date" => String),
    "constants" => ("num_days" => Integer, "dt" => Integer, "dx" => Real, "dy" => Real, "num_parc" => Integer, "atm_mass" => Real),
    "paths" => ("path_data1" => String, "path_data2" => String, "path_mask" => String, "path_obs" => String),
    "names" => ("obs_evap" => String, "obs_heat" => String, "obs_prec" => String),
    "flags" => ("forward" => Bool, "track_heat" => Bool, "obs_era5" => Bool, "save_q" => Bool, "save_diagnose" => Bool),
)

"""
    read_config(path)

Settings from the TOML namelist at `path`. Every key of `NAMELIST` must be there with a value of its type, no `set_...`
placeholder may be left, dates are `yyyy-mm-dd_HH` and `obs_era5` must be true. All problems are reported in one error.
"""
function read_config(path)
    config = TOML.parsefile(path)
    problems = String[]
    for (section, keys) in NAMELIST, (key, type) in keys
        value = get(get(config, section, Dict()), key, nothing)
        if isnothing(value)
            push!(problems, "[$section] $key is missing")
        elseif value isa String && occursin(r"(^|/)set_", value)
            push!(problems, "[$section] $key is still the placeholder \"$value\"")
        elseif !(value isa type) || (value isa Bool && type !== Bool)
            push!(problems, "[$section] $key = $(repr(value)) is not of type $type")
        elseif section == "dates" && isnothing(tryparse(DateTime, value, DATEFORMAT))
            push!(problems, "[$section] $key = \"$value\" is not a date of the form yyyy-mm-dd_HH")
        end
    end
    get(get(config, "flags", Dict()), "obs_era5", true) === false && push!(problems, "[flags] obs_era5 = false is not supported")
    isempty(problems) || error("namelist $path:\n  " * join(problems, "\n  "))
    return config
end

"Longitude and latitude `(lon, lat)` of every sink cell (mask equal to 1) on a 0.25 degree grid."
function read_mask(path)
    lon, lat, mask = NCDataset(path) do ds
        m = ds["mask"]
        if ndims(m) == 3
            mask = m[:, :, 1]
        else
            @info "Source mask has no temporal dimension"
            mask = m[:, :]
        end
        ds["lon"][:], ds["lat"][:], mask
    end
    lons = -180:0.25:179.75
    lats = -90:0.25:90
    nearest(grid, x) = argmin(abs.(grid .- x))
    ix = [nearest(lon, x) for x in lons]
    iy = [nearest(lat, y) for y in lats]
    return [(lons[i], lats[j]) for j in eachindex(lats) for i in eachindex(lons) if isequal(mask[ix[i], iy[j]], 1)]
end

"Paths of the partoutput files of the given dates."
partoutput_files(path, dates) =
    [joinpath(path, "partoutputs", "partoutput_" * stamp(d) * ".nc") for d in dates]

"""
    read_window(files, parcels, names)

Read the variables `names` of the parcels `parcels` (an index range) from every file. Returns a named tuple of
`(parcels, files)` arrays. Longitudes are wrapped to -180 to 180 and NaN latitudes are set to -1.
"""
function read_window(files, parcels, names)
    data = [zeros(Float32, length(parcels), length(files)) for _ in names]
    for (i, file) in enumerate(files)
        NCDataset(file) do ds
            for (a, name) in zip(data, names)
                a[:, i] = ds[name].var[1, parcels]
            end
        end
    end
    window = NamedTuple{Symbol.(Tuple(names))}(data)
    haskey(window, :lat) && replace!(window.lat, NaN32 => -1.0f0)
    haskey(window, :lon) && (@. window.lon = ifelse(window.lon > 180, window.lon - 360, window.lon))
    return window
end

"""
    read_obs(file, name, time, g)

Field `name` of an ERA5 file at `time`, shifted to longitudes -180 to 180 and interpolated linearly onto the output
grid `g` (see `grid_definition`). Returns an `(nx, ny)` array.
"""
function read_obs(file, name, time, g)
    NCDataset(file) do ds
        it = findfirst(==(time), ds["time"][:])
        isnothing(it) && error("time $time is not in $file")
        lon = mod.(ds["longitude"][:] .+ 180, 360) .- 180
        lat = ds["latitude"][:]
        field = nomissing(ds[name][:, :, it], NaN)
        ilon, ilat = sortperm(lon), sortperm(lat)
        itp = linear_interpolation((lon[ilon], lat[ilat]), field[ilon, ilat])
        return [itp(x, y) for x in g.lons, y in g.lats]
    end
end

"""
    read_observations(config, time)

Observed evaporation, precipitation and (with `track_heat`) sensible heat at `time`, as `(nx, ny)` Float32 arrays
`(e_obs, p_obs, h_obs)`. Evaporation and heat are >= 0 and precipitation <= 0; `h_obs` is `nothing` without heat.
"""
function read_observations(config, time)
    path, names = config["paths"]["path_obs"], config["names"]
    dt, g = config["constants"]["dt"], grid_definition(config)
    obs(file, name) = Float32.(read_obs(joinpath(path, file * ".nc"), name, time, g))
    e_obs = max.(obs(names["obs_evap"], "e") .* -1000.0f0, 0.0f0)
    p_obs = min.(obs(names["obs_prec"], "tp") .* -1000.0f0, 0.0f0)
    h_obs = config["flags"]["track_heat"] ? max.(obs(names["obs_heat"], "sshf") ./ Float32(-dt * 3600), 0.0f0) : nothing
    return e_obs, p_obs, h_obs
end

"Define the dimensions `lon`, `lat`, `time` (one step, `time`) and their coordinates of grid `g` in the dataset `ds`."
function define_grid!(ds, g, time)
    defDim(ds, "lon", g.nx)
    defDim(ds, "lat", g.ny)
    defDim(ds, "time", 1)
    defVar(ds, "lon", Float32.(g.lons), ("lon",))
    defVar(ds, "lat", Float32.(g.lats), ("lat",))
    defVar(ds, "time", [time], ("time",))
end

"""
    write_diagnosis(file, config, time, d)

Write the gridded fields of `diagnose` (`d`) for `time` to a NetCDF file with dimensions `(lon, lat, time)`: `E`, `P`,
`IWV` and, if present, `H`, each with the number of parcels per cell in `*_n_part` (`n_part` for `IWV`).
"""
function write_diagnosis(file, config, time, d)
    g = grid_definition(config)
    fields = (E=d.evap, E_n_part=d.n_evap, P=d.prec, P_n_part=d.n_prec, H=d.heat, H_n_part=d.n_heat, IWV=d.iwv, n_part=d.n_part)
    mkpath(dirname(file))
    NCDataset(file, "c") do ds
        define_grid!(ds, g, time)
        for (name, a) in pairs(fields)
            isnothing(a) || (defVar(ds, String(name), eltype(a), ("lon", "lat", "time"))[:, :, 1] = a)
        end
    end
end

"""
    write_corrected(file, raw_file, parcels, extras)

Write the file that the tracking reads for one time step: `lon`, `lat` and `sh` of the parcels `parcels` copied from
`raw_file`, and the variables of the named tuple `extras` (`Delta_q`, `theta`, `Delta_theta`; entries that are
`nothing` are skipped). Same layout as the raw files, so `read_window` reads it.
"""
function write_corrected(file, raw_file, parcels, extras)
    mkpath(dirname(file))
    NCDataset(raw_file) do raw
        NCDataset(file, "c") do ds
            defDim(ds, "time", 1)
            defDim(ds, "particle", length(parcels))
            for name in ("lon", "lat", "sh")
                defVar(ds, name, Float32, ("time", "particle"))[1, :] = raw[name].var[1, parcels]
            end
            for (name, a) in pairs(extras)
                isnothing(a) || (defVar(ds, String(name), Float32, ("time", "particle"))[1, :] = a)
            end
        end
    end
end

"""
    write_attribution(file, config, time, tlevel, name, field)

Write one attribution field `(nx, ny, days)` (`name` is `E2P`, `E2Q` or `H2T`) for `time` to a NetCDF file with
dimensions `(lon, lat, tlevel, time)`; `tlevel` are the days the contributions come from.
"""
function write_attribution(file, config, time, tlevel, name, field)
    g = grid_definition(config)
    mkpath(dirname(file))
    NCDataset(file, "c") do ds
        define_grid!(ds, g, time)
        defDim(ds, "tlevel", length(tlevel))
        defVar(ds, "tlevel", tlevel, ("tlevel",))
        defVar(ds, name, Float32, ("lon", "lat", "tlevel", "time"))[:, :, :, 1] = field
    end
end
