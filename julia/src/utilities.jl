"Midpoint of two positions on the sphere, in degrees. Returns `(lon, lat)`."
function midpoints(lona, lonb, lata, latb)
    xmidp = @. 0.5f0 * cosd(latb) * cosd(lonb)
    xmidp = @. xmidp + 0.5f0 * cosd(lata) * cosd(lona)
    ymidp = @. 0.5f0 * cosd(latb) * sind(lonb)
    ymidp = @. ymidp + 0.5f0 * cosd(lata) * sind(lona)
    zmidp = @. 0.5f0 * (sind(latb) + sind(lata))
    lonmidp = @. atand(ymidp, xmidp)
    latmidp = @. atand(zmidp, sqrt(xmidp * xmidp + ymidp * ymidp))
    return lonmidp, latmidp
end

"Potential temperature from specific humidity, temperature (K) and density."
function only_theta(qv::AbstractArray{T}, temp, rho) where {T<:AbstractFloat}
    R_d, R_w = T(287.057), T(461.5)
    ε = R_d / R_w
    w = @. qv / (1 - qv)
    tv = @. temp / (1 - w * (1 - ε) / (w + ε))
    p_hpa = @. rho * R_d * tv / 100
    return @. temp * (1000 / p_hpa)^(T(0.2854) * (1 - T(0.00028) * w))
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
#     theta = @. temp * (1000 / p_hpa)^(0.2854 * (1 - 0.00028 * w))
#     return RH, theta
# end

"""
    ongrid(nx, ny, x_index, y_index, f, mass)

Add `f * mass` of each parcel to its grid cell, for every time step. `f`, `x_index` and
`y_index` are `(parcels, times)`; returns a Float32 array `(nx, ny, times)`.
"""
function ongrid(nx, ny, x_index, y_index, f, mass)
    numpart, ntimes = size(f)
    ms = zeros(Float32, nx, ny, ntimes)
    for t in 1:ntimes
        for j in 1:numpart
            # ms[x_index[j, t], y_index[j, t], t] = ms[x_index[j, t], y_index[j, t], t] + f[j, t] * mass[j, t]
            ms[x_index[j, t], y_index[j, t], t] = ms[x_index[j, t], y_index[j, t], t] + f[j, t] * mass
        end
    end
    return ms
end

"""
    ongrid1d(nx, ny, x_index, y_index, f, mass)

One time step of `ongrid`, for all parcels but the last. Returns the summed `f * mass` and
the number of parcels per cell, both `(nx, ny)`.
"""
function ongrid1d(nx, ny, x_index, y_index, f, mass)
    var1 = zeros(Float32, nx, ny)
    var2 = zeros(Float32, nx, ny)
    for j in 1:length(f)-1
        # var1[x_index[j], y_index[j]] = var1[x_index[j], y_index[j]] + f[j] * mass[j]
        var1[x_index[j], y_index[j]] = var1[x_index[j], y_index[j]] + f[j] * mass
        var2[x_index[j], y_index[j]] = var2[x_index[j], y_index[j]] + 1
    end
    return var1, var2
end

"""
    hmin(nx, ny, x_index, y_index, h1, h2)

Lowest `h1` and `h2` per grid cell for all parcels but the last, 30000 where a cell has none.
"""
function hmin(nx, ny, x_index, y_index, h1, h2)
    var1 = fill(30000.0f0, nx, ny)
    var2 = fill(30000.0f0, nx, ny)
    for j in 1:length(h1)-1
        var1[x_index[j], y_index[j]] = min(var1[x_index[j], y_index[j]], h1[j])
        var2[x_index[j], y_index[j]] = min(var2[x_index[j], y_index[j]], h2[j])
    end
    return var1, var2
end
