"Midpoint of two positions on the sphere, in degrees. Returns `(lon, lat)`."
function midpoint(lona, lonb, lata, latb)
    rad(x) = x * Float32(π) / 180.0f0
    deg(x) = x * 180.0f0 / Float32(π)
    xmidp = 0.5f0 * cos(rad(latb)) * cos(rad(lonb))
    xmidp = xmidp + 0.5f0 * cos(rad(lata)) * cos(rad(lona))
    ymidp = 0.5f0 * cos(rad(latb)) * sin(rad(lonb))
    ymidp = ymidp + 0.5f0 * cos(rad(lata)) * sin(rad(lona))
    zmidp = 0.5f0 * (sin(rad(latb)) + sin(rad(lata)))
    return deg(atan(ymidp, xmidp)), deg(atan(zmidp, sqrt(xmidp * xmidp + ymidp * ymidp)))
end

"`midpoint` of every pair of positions in the arrays. Returns `(lon, lat)` arrays."
function midpoints(lona, lonb, lata, latb)
    m = midpoint.(lona, lonb, lata, latb)
    return first.(m), last.(m)
end

"""
    cells(g, lona, lonb, lata, latb)

Grid cells `(x, y)` (Int32) on grid `g` of the midpoints of the positions `a` and `b`; the column at 180 degrees wraps
to the first. NaN positions get `typemin(Int32)`, outside the grid.
"""
function cells(g, lona, lonb, lata, latb)
    x, y = similar(lona, Int32), similar(lona, Int32)
    dx, dy = Float32(g.dx), Float32(g.dy)
    for i in eachindex(x, lona, lonb, lata, latb)
        lon, lat = midpoint(lona[i], lonb[i], lata[i], latb[i])
        xi = gridindex(lon + 180, dx)
        x[i] = xi == g.nx + 1 ? 1 : xi
        y[i] = gridindex(lat + 90, dy)
    end
    return x, y
end

"Potential temperature from specific humidity, temperature (K) and density."
function only_theta(qv::AbstractArray{T}, temp, rho) where {T<:AbstractFloat}
    R_d, R_w = T(287.057), T(461.5)
    ε = R_d / R_w
    w = @. qv / (1 - qv)
    tv = @. temp / (1 - w * (1 - ε) / (w + ε))
    p_hpa = @. rho * R_d * tv / 100
    return @. temp * (1000 / p_hpa)^(T(0.2854) * (1 - T(0.28) * w))
end

"""
    ongrid(nx, ny, x_index, y_index, f, mass)

Add `f * mass` of each parcel to its grid cell, for every time step. Parcels outside the grid are skipped. `f`, `x_index` and
`y_index` are `(parcels, times)`; returns a Float32 array `(nx, ny, times)`.
"""
function ongrid(nx, ny, x_index, y_index, f, mass)
    numpart, ntimes = size(f)
    ms = zeros(Float32, nx, ny, ntimes)
    for t in 1:ntimes
        for j in 1:numpart
            # ms[x_index[j, t], y_index[j, t], t] = ms[x_index[j, t], y_index[j, t], t] + f[j, t] * mass[j, t]
            if 1 <= x_index[j, t] <= nx && 1 <= y_index[j, t] <= ny
                ms[x_index[j, t], y_index[j, t], t] = ms[x_index[j, t], y_index[j, t], t] + f[j, t] * mass
            end
        end
    end
    return ms
end

"""
    ongrid1d(nx, ny, x_index, y_index, f, mass)

One time step of `ongrid`. Returns the summed `f * mass` and
the number of parcels per cell, both `(nx, ny)`.
"""
function ongrid1d(nx, ny, x_index, y_index, f, mass)
    var1 = zeros(Float32, nx, ny)
    var2 = zeros(Float32, nx, ny)
    for j in eachindex(f)
        # var1[x_index[j], y_index[j]] = var1[x_index[j], y_index[j]] + f[j] * mass[j]
        if 1 <= x_index[j] <= nx && 1 <= y_index[j] <= ny
            var1[x_index[j], y_index[j]] = var1[x_index[j], y_index[j]] + f[j] * mass
            var2[x_index[j], y_index[j]] = var2[x_index[j], y_index[j]] + 1
        end
    end
    return var1, var2
end

"""
    grid_definition(config)

Output grid from the settings: size `nx`, `ny`, spacing `dx`, `dy`, cell coordinates `lons`, `lats` (from -180 and -90),
mass per parcel and cell areas `(nx, ny)` in m^2 (NaN where the area is 0).
"""
function grid_definition(config)
    c = config["constants"]
    dx, dy = c["dx"], c["dy"]
    nx, ny = trunc(Int, 360 / dx), trunc(Int, 180 / dy) + 1
    lons, lats = range(-180, step=dx, length=nx), range(-90, step=dy, length=ny)
    lat = Float32.(lats)'
    lats1 = @. (lat + Float32(dy / 2)) * Float32(π) / 180.0f0
    lats2 = @. (lat - Float32(dy / 2)) * Float32(π) / 180.0f0
    area = Float32(dx * π / 180 * 6371000.0^2) .* (sin.(lats1) .- sin.(lats2))
    areas = repeat(ifelse.(area .== 0, NaN32, area), nx)
    return (; nx, ny, dx, dy, lons, lats, mass=c["atm_mass"] / c["num_parc"], areas)
end

"Grid column or row (1-based, Int32) of `x`, a coordinate counted from the grid start, for spacing `d`. NaN falls outside the grid."
gridindex(x, d) = isfinite(x) ? round(Int32, x / d) + Int32(1) : typemin(Int32)

"Replace NaN by 0."
nan0(a) = replace(a, NaN32 => 0.0f0)

"Gridded water mass (kg per cell) as a depth in mm on grid `g`; 0 where the cell has no area."
to_mm(g, kg) = nan0(kg .* Float32(1000 / 997) ./ g.areas)

"Gridded heat (theta change times mass per cell) as a flux in W/m^2 over `seconds` on grid `g`; 0 where the cell has no area."
to_wm2(g, heat, seconds) = nan0(heat .* Float32(1005.7) ./ Float32(seconds) ./ g.areas)
