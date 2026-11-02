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
