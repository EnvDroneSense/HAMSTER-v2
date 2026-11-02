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
    [joinpath(path, "partoutputs", "partoutput_" * Dates.format(d, "yyyymmddHHMMSS") * ".nc") for d in dates]

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
    read_obs(file, name, time, dx, dy)

Field `name` of an ERA5 file at `time`, shifted to longitudes -180 to 180 and interpolated linearly onto the output
grid. Returns an `(nx, ny)` array.
"""
function read_obs(file, name, time, dx, dy)
    lons = range(-180, step=dx, length=round(Int, 360 / dx))
    lats = range(-90, step=dy, length=round(Int, 180 / dy) + 1)
    NCDataset(file) do ds
        it = findfirst(==(time), ds["time"][:])
        isnothing(it) && error("time $time is not in $file")
        lon = mod.(ds["longitude"][:] .+ 180, 360) .- 180
        lat = ds["latitude"][:]
        field = nomissing(ds[name][:, :, it], NaN)
        ilon, ilat = sortperm(lon), sortperm(lat)
        itp = linear_interpolation((lon[ilon], lat[ilat]), field[ilon, ilat])
        return [itp(x, y) for x in lons, y in lats]
    end
end

"""
    read_observations(config, time)

Observed evaporation, precipitation and (with `track_heat`) sensible heat at `time`, as `(nx, ny)` Float32 arrays
`(e_obs, p_obs, h_obs)`. Evaporation and heat are >= 0 and precipitation <= 0; `h_obs` is `nothing` without heat.
"""
function read_observations(config, time)
    path, names = config["paths"]["path_obs"], config["names"]
    dt, dx, dy = config["constants"]["dt"], config["constants"]["dx"], config["constants"]["dy"]
    obs(file, name) = Float32.(read_obs(joinpath(path, file * ".nc"), name, time, dx, dy))
    e_obs = max.(obs(names["obs_evap"], "e") .* -1000.0f0, 0.0f0)
    p_obs = min.(obs(names["obs_prec"], "tp") .* -1000.0f0, 0.0f0)
    h_obs = config["flags"]["track_heat"] ? max.(obs(names["obs_heat"], "sshf") ./ Float32(-dt * 3600), 0.0f0) : nothing
    return e_obs, p_obs, h_obs
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
        defDim(ds, "lon", g.nx)
        defDim(ds, "lat", g.ny)
        defDim(ds, "time", 1)
        defVar(ds, "lon", Float32.(range(-180, step=g.dx, length=g.nx)), ("lon",))
        defVar(ds, "lat", Float32.(range(-90, step=g.dy, length=g.ny)), ("lat",))
        defVar(ds, "time", [time], ("time",))
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
        defDim(ds, "lon", g.nx)
        defDim(ds, "lat", g.ny)
        defDim(ds, "tlevel", length(tlevel))
        defDim(ds, "time", 1)
        defVar(ds, "lon", Float32.(range(-180, step=g.dx, length=g.nx)), ("lon",))
        defVar(ds, "lat", Float32.(range(-90, step=g.dy, length=g.ny)), ("lat",))
        defVar(ds, "tlevel", tlevel, ("tlevel",))
        defVar(ds, "time", [time], ("time",))
        defVar(ds, name, Float32, ("lon", "lat", "tlevel", "time"))[:, :, :, 1] = field
    end
end
