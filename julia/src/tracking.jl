"Set a value to 0 if it is negative or not finite."
clean(x) = (x < 0 || !isfinite(x)) ? zero(x) : x

"""
    adjust_losses(q, dq)

Contribution of each time interval's moisture (or heat) gain to the parcel at the end of the trajectory. `q` and
`dq` are `(parcels, intervals)`: the quantity at the start of each interval and its change. A gain is kept as
is; a loss reduces the earlier contributions in proportion to what the parcel loses. The last interval is left 0.
Returns Float64.
"""
function adjust_losses(q, dq)
    n, m = size(dq)
    out = zeros(Float64, n, m)
    for j in 1:n, t in 1:m-1
        d = dq[j, t]
        if d > 0
            out[j, t] = d
        elseif d < 0
            r = (q[j, t] + d) / q[j, t]
            for s in 1:t-1
                out[j, s] *= r
            end
        end
        for s in 1:t
            out[j, s] = clean(out[j, s])
        end
    end
    return out
end

"""
    arrival_precipitation(e2q, dq, q)

Contribution to the precipitation in the last interval for the parcels that lose moisture there: `e2q` is the
`adjust_losses` result of these parcels, `dq` and `q` their change and moisture in the last interval. The
contributions are scaled so that they add up to the parcel's whole loss.
"""
function arrival_precipitation(e2q, dq, q)
    e = @. -e2q * (dq / q)
    fraction = Float64.(dq) ./ vec(sum(e, dims=2))
    return clean.(@. -e * fraction)
end

"""
    on_sink(sink, lon, lat)

Which parcels (rows of `lon` and `lat`, `(parcels, times)`) are in the sink between the last two time steps: their
midpoint, rounded to 0.25 degrees, is in the set `sink` of `(lon, lat)` tuples (see `read_mask`).
"""
function on_sink(sink, lon, lat)
    lonmid, latmid = midpoints(lon[:, end], lon[:, end-1], lat[:, end], lat[:, end-1])
    snap(x) = Float64(round(x * 4.0f0) / 4.0f0) + 0.0   # + 0.0 turns -0.0 into 0.0 so that it matches 0.0 in the sink set
    return [(snap(x), snap(y)) in sink for (x, y) in zip(lonmid, latmid)]
end

"""
    attribute(config, sink, datesmid, lon, lat, qv, Delta_q, theta, Delta_theta)

Moisture (and heat) sources of the parcels that are in the `sink` at the end of the window. The arrays are
`(parcels, times)`; `datesmid` are the dates in the middle of the intervals. Returns the named tuple `(E2P, E2Q, H2T)` of
daily `(nx, ny, days)` fields in mm (W/m^2 for `H2T`); `E2Q` and `H2T` are `nothing` when the flags are off.
"""
function attribute(config, sink, datesmid, lon, lat, qv, Delta_q, theta, Delta_theta)
    g, flags = grid_definition(config), config["flags"]
    sel = on_sink(sink, lon, lat)
    q, dq = qv[sel, :], Delta_q[sel, 2:end]
    m = size(dq, 2)

    e2q = adjust_losses(q[:, 1:m], dq)
    arrival = dq[:, m] .< 0
    e2p = arrival_precipitation(e2q[arrival, :], dq[arrival, m], q[arrival, m])

    lonmid, latmid = midpoints(lon[sel, 2:end], lon[sel, 1:end-1], lat[sel, 2:end], lat[sel, 1:end-1])
    x = gridindex.(lonmid .+ 180, Float32(g.dx))
    x = @. ifelse(x == g.nx + 1, 1, x)
    y = gridindex.(latmid .+ 90, Float32(g.dy))

    E2P = ongrid(g.nx, g.ny, x[arrival, :], y[arrival, :], e2p, g.mass)
    E2Q = flags["save_e2q"] ? ongrid(g.nx, g.ny, x, y, e2q, g.mass) : nothing
    H2T = nothing
    if flags["track_heat"]
        h2t = adjust_losses(theta[sel, 1:m], Delta_theta[sel, 2:end])
        H2T = ongrid(g.nx, g.ny, x, y, h2t, g.mass)
    end

    days = Date.(datesmid)
    daily(a) = stack(d -> dropdims(sum(a[:, :, days .== d], dims=3), dims=3), unique(days))
    mm(a) = nan0(daily(a) .* Float32(1000 / 997) ./ g.areas)
    wm2(a) = nan0(daily(a) .* Float32(1005.7) ./ Float32(24 * 3600) ./ g.areas)
    return (; E2P=mm(E2P), E2Q=isnothing(E2Q) ? nothing : mm(E2Q), H2T=isnothing(H2T) ? nothing : wm2(H2T))
end

"Parcel ranges of `n` blocks of equal size; the last block takes the remainder."
function parcel_blocks(num_parc, n)
    len = num_parc ÷ n
    return [(i-1)*len+1:(i == n ? num_parc : i * len) for i in 1:n]
end

"Drop the first time step of `window` (a named tuple of `(parcels, times)` arrays) and append the new one."
shift_window(window, new) = map((a, b) -> hcat(a[:, 2:end], b), window, new)

"Sum of the fields of the block results `rs` (the gather and sum over the ranks); `nothing` where a field is off."
combine(rs) = map((fields...) -> isnothing(first(fields)) ? nothing : sum(fields), rs...)

"""
    track(config; blocks=1)

Attribution for every arrival date from `start_date + dt` to `end_date`: for each date, `attribute` is run on the
`num_days` before it and `attribution_{e2p,e2q,h2t}_<arrival>.nc` is written to `path_data2/attribution`. The parcels
are split into `blocks` parts that are processed one after the other and summed. The first date reads all files;
every later date reads only its newest file and shifts the window.
"""
function track(config; blocks=1)
    c, d, flags, paths = config["constants"], config["dates"], config["flags"], config["paths"]
    dt, half = c["dt"], Second(c["dt"] * 1800)
    format = dateformat"yyyy-mm-dd_HH"
    sink = Set(read_mask(joinpath(paths["path_mask"], "mask.nc")))
    names = ["lon", "lat", "sh", "Delta_q"]
    flags["track_heat"] && append!(names, ["theta", "Delta_theta"])
    parcels = parcel_blocks(c["num_parc"], blocks)

    window = nothing
    for arrival in DateTime(d["start_date"], format)+Hour(dt):Hour(dt):DateTime(d["end_date"], format)
        starting = arrival - Day(c["num_days"])
        dates = isnothing(window) ? (starting:Hour(dt):arrival) : [arrival]
        datesmid = starting+half:Hour(dt):arrival-half
        files = partoutput_files(joinpath(paths["path_data1"], "corrected"), dates)
        new = [read_window(files, p, names) for p in parcels]
        window = isnothing(window) ? new : shift_window.(window, new)

        r = combine([attribute(config, sink, datesmid, w.lon, w.lat, w.sh, w.Delta_q, get(w, :theta, nothing), get(w, :Delta_theta, nothing)) for w in window])
        tlevel = DateTime.(unique(Date.(datesmid)))
        stamp = Dates.format(arrival, "yyyymmddHHMMSS")
        for (name, field) in pairs(r)
            isnothing(field) || write_attribution(joinpath(paths["path_data2"], "attribution", "attribution_$(lowercase(String(name)))_$stamp.nc"), config, arrival - half, tlevel, String(name), field)
        end
    end
end
