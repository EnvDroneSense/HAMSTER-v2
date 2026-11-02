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
