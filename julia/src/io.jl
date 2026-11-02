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

"""
    read_window(files, parcels, track_heat)

Read the parcels `parcels` (an index range) from every file. Returns `(lon, lat, qv, Delta_q, theta, Delta_theta)`,
each `(parcels, files)`; `theta` and `Delta_theta` are `nothing` unless `track_heat`. Longitudes are in -180 to 180
and NaN latitudes are set to -1.
"""
function read_window(files, parcels, track_heat)
    names = track_heat ? ("lon", "lat", "sh", "Delta_q", "theta", "Delta_theta") : ("lon", "lat", "sh", "Delta_q")
    data = [zeros(Float32, length(parcels), length(files)) for _ in names]
    for (i, file) in enumerate(files)
        NCDataset(file) do ds
            for (a, name) in zip(data, names)
                a[:, i] = ds[name].var[parcels, 1]
            end
        end
    end
    lon, lat, qv, Delta_q = data
    replace!(lat, NaN32 => -1.0f0)
    @. lon = ifelse(lon > 180, lon - 360, lon)
    theta, Delta_theta = track_heat ? (data[5], data[6]) : (nothing, nothing)
    return lon, lat, qv, Delta_q, theta, Delta_theta
end
