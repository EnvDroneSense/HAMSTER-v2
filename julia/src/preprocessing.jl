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

"""
    preprocess_window(config, arrival, first)

Bias-correct the changes between `arrival - dt` and `arrival`: read the two partoutputs, correct the changes against
the observations and write the corrected files to `path_data1/corrected` (and the diagnosis to `path_data1/diagnosis`
with `save_diagnose`). The raw files are not changed. `first` also writes the file of `arrival - dt`, with zero changes.
"""
function preprocess_window(config, arrival, first)
    c, flags, path = config["constants"], config["flags"], config["paths"]["path_data1"]
    dt, parcels, heat = c["dt"], 1:c["num_parc"], config["flags"]["track_heat"]
    before, middle = arrival - Hour(dt), arrival - Second(dt * 1800)
    raw = partoutput_files(path, [before, arrival])
    out = partoutput_files(joinpath(path, "corrected"), [before, arrival])

    w = read_window(raw, parcels, heat ? ["lon", "lat", "sh", "T", "rho"] : ["lon", "lat", "sh"])
    theta = heat ? only_theta(w.sh, w.T, w.rho) : nothing
    Delta_q = delta(w.sh)
    Delta_theta = heat ? delta(theta) : nothing

    d = diagnose(config, w.lon, w.lat, w.sh, Delta_q, Delta_theta)
    flags["save_diagnose"] && write_diagnosis(joinpath(path, "diagnosis", "diagnosis_" * Dates.format(arrival, "yyyymmddHHMMSS") * ".nc"), config, middle, d)
    Delta_q, Delta_theta = correct_deltas(d, Delta_q, Delta_theta, read_observations(config, middle)...)

    column(a, i) = isnothing(a) ? nothing : a[:, i]
    if first
        zero_q = zeros(Float32, length(parcels))
        write_corrected(out[1], raw[1], parcels, (; Delta_q=zero_q, theta=column(theta, 1), Delta_theta=heat ? zero_q : nothing))
    end
    write_corrected(out[2], raw[2], parcels, (; Delta_q, theta=column(theta, 2), Delta_theta))
end

"""
Bias-correct all time windows (see `preprocess_window`): arrival dates from `start_date - num_days + dt` to `end_date`;
forward from `start_date + dt` to `end_date + num_days`.
"""
function preprocess(config)
    c, d = config["constants"], config["dates"]
    format = dateformat"yyyy-mm-dd_HH"
    forward = config["flags"]["forward"]
    start, stop = DateTime(d["start_date"], format), DateTime(d["end_date"], format)
    first_date = forward ? start + Hour(c["dt"]) : start - Day(c["num_days"]) + Hour(c["dt"])
    for arrival in first_date:Hour(c["dt"]):(forward ? stop + Day(c["num_days"]) : stop)
        preprocess_window(config, arrival, arrival == first_date)
    end
end
