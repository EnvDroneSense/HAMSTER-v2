"Midpoint of two positions on the sphere, in degrees. Returns `(lon, lat)`."
function midpoints(lona, lonb, lata, latb)
    rad(x) = x * Float32(π) / 180.0f0
    deg(x) = x * 180.0f0 / Float32(π)
    xmidp = @. 0.5f0 * cos(rad(latb)) * cos(rad(lonb))
    xmidp = @. xmidp + 0.5f0 * cos(rad(lata)) * cos(rad(lona))
    ymidp = @. 0.5f0 * cos(rad(latb)) * sin(rad(lonb))
    ymidp = @. ymidp + 0.5f0 * cos(rad(lata)) * sin(rad(lona))
    zmidp = @. 0.5f0 * (sin(rad(latb)) + sin(rad(lata)))
    lonmidp = @. deg(atan(ymidp, xmidp))
    latmidp = @. deg(atan(zmidp, sqrt(xmidp * xmidp + ymidp * ymidp)))
    return lonmidp, latmidp
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

# function RH_theta(qv, temp, rho)
#     R_d, R_w = 287.057, 461.5
#     ε = R_d / R_w
#     w = @. qv / (1 - qv)
#     tv = @. temp / (1 - w * (1 - ε) / (w + ε))
#     p_pa = @. rho * R_d * tv
#     e = @. p_pa * w / (w + ε)
#     e_s = @. 611.2 * exp(17.67 * (temp - 273.15) / (temp - 273.15 + 243.5))
#     RH = @. 100 * e / e_s
#     p_hpa = @. p_pa / 100
#     theta = @. temp * (1000 / p_hpa)^(0.2854 * (1 - 0.28 * w))
#     return RH, theta
# end

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
    hmin(nx, ny, x_index, y_index, h1, h2)

Lowest `h1` and `h2` per grid cell, 30000 where a cell has none. Parcels outside the grid are skipped.
"""
function hmin(nx, ny, x_index, y_index, h1, h2)
    var1 = fill(30000.0f0, nx, ny)
    var2 = fill(30000.0f0, nx, ny)
    for j in eachindex(h1)
        if 1 <= x_index[j] <= nx && 1 <= y_index[j] <= ny
            var1[x_index[j], y_index[j]] = min(var1[x_index[j], y_index[j]], h1[j])
            var2[x_index[j], y_index[j]] = min(var2[x_index[j], y_index[j]], h2[j])
        end
    end
    return var1, var2
end

"Grid size, spacing, mass per parcel and cell areas `(nx, ny)` in m^2 (NaN where the area is 0) from the settings."
function grid_definition(config)
    c = config["constants"]
    dx, dy = c["dx"], c["dy"]
    nx, ny = trunc(Int, 360 / dx), trunc(Int, 180 / dy) + 1
    lat = Float32.(range(-90, step=dy, length=ny))'
    lats1 = @. (lat + Float32(dy / 2)) * Float32(π) / 180.0f0
    lats2 = @. (lat - Float32(dy / 2)) * Float32(π) / 180.0f0
    area = Float32(dx * π / 180 * 6371000.0^2) .* (sin.(lats1) .- sin.(lats2))
    areas = repeat(ifelse.(area .== 0, NaN32, area), nx)
    return (; nx, ny, dx, dy, mass=c["atm_mass"] / c["num_parc"], areas)
end

"Grid column or row (1-based) of `x`, a coordinate counted from the grid start, for spacing `d`. NaN falls outside the grid."
gridindex(x, d) = isfinite(x) ? round(Int, x / d) + 1 : typemin(Int)

"Replace NaN by 0."
nan0(a) = replace(a, NaN32 => 0.0f0)
