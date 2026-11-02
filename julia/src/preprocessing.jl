"Change between the two time steps (columns) of `x`; changes of at most 1e-4 in size are set to 0."
function delta(x)
    d = x[:, 2] .- x[:, 1]
    return @. ifelse(abs(d) <= 1.0f-4, zero(d), d)
end

"""
    diagnose(config, lon, lat, qv, Delta_q, Delta_theta)

Grid the evaporation, precipitation and heat that the parcels diagnose between the two time steps (columns) of
`lon`, `lat` and `qv`. Returns the parcels' grid indices and the `(nx, ny)` fields `iwv`, `evap`, `prec`, `heat`
(mm, mm, mm and W/m^2; `heat` is `nothing` without `Delta_theta`) with the number of parcels in each cell (`n_*`).
"""
function diagnose(config, lon, lat, qv, Delta_q, Delta_theta)
    g = grid_definition(config)
    lonmid, latmid = midpoints(lon[:, 2], lon[:, 1], lat[:, 2], lat[:, 1])
    x = gridindex.(lonmid .+ 180, Float32(g.dx))
    x = @. ifelse(x == g.nx + 1, 1, x)
    y = gridindex.(latmid .+ 90, Float32(g.dy))

    grid(f, keep) = ongrid1d(g.nx, g.ny, x[keep], y[keep], f[keep], g.mass)
    nan0(a) = replace(a, NaN32 => 0.0f0)
    iwv, n_part = ongrid1d(g.nx, g.ny, x, y, vec(mean(qv, dims=2)), g.mass)
    evap, n_evap = grid(Delta_q, Delta_q .> 0)
    prec, n_prec = grid(Delta_q, Delta_q .< 0)
    mm(kg) = nan0(kg .* Float32(1000 / 997) ./ g.areas)
    heat = n_heat = nothing
    if !isnothing(Delta_theta)
        heat, n_heat = grid(Delta_theta, Delta_theta .> 0)
        heat = nan0(heat .* Float32(1005.7) ./ Float32(config["constants"]["dt"] * 3600) ./ g.areas)
    end
    return (; x, y, iwv=mm(iwv), n_part, evap=mm(evap), n_evap, prec=mm(prec), n_prec, heat, n_heat)
end

"The values of `field[x, y]` for each parcel; NaN for parcels outside the grid."
sample(field, x, y) = [1 <= i <= size(field, 1) && 1 <= j <= size(field, 2) ? field[i, j] : NaN32 for (i, j) in zip(x, y)]

"""
    correct(delta, up, down)

Scale positive `delta` by `up` and negative `delta` by `down`; zero stays zero and non-finite results become 0.
"""
function correct(delta, up, down)
    new = @. ifelse(delta > 0, delta * up, ifelse(delta < 0, delta * down, zero(delta)))
    return @. ifelse(isfinite(new), new, zero(new))
end

"""
    correct_deltas(d, Delta_q, Delta_theta, e_obs, p_obs, h_obs)

Bias-correct the humidity and heat changes: scale each parcel's change by the ratio of observed to diagnosed
evaporation, precipitation or heat in its grid cell (`d` is the result of `diagnose`). Negative theta changes are kept.
"""
function correct_deltas(d, Delta_q, Delta_theta, e_obs, p_obs, h_obs)
    ratio(obs, diag) = sample(obs, d.x, d.y) ./ sample(diag, d.x, d.y)
    Delta_q = correct(Delta_q, ratio(e_obs, d.evap), ratio(p_obs, d.prec))
    Delta_theta = isnothing(Delta_theta) ? nothing : correct(Delta_theta, ratio(h_obs, d.heat), 1)
    return Delta_q, Delta_theta
end
