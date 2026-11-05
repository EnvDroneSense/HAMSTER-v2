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
    out = zeros(Float64, m, n)   # intervals first, so that each parcel's contributions are contiguous
    @inbounds for j in 1:n, t in 1:m-1
        d = dq[j, t]
        if d > 0
            out[t, j] = clean(d)
        elseif d < 0
            r = (q[j, t] + d) / q[j, t]
            for s in 1:t-1
                out[s, j] = clean(out[s, j] * r)
            end
        end
    end
    return permutedims(out)
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
    follow_losses(q, dq)

Forward counterpart of `adjust_losses`: the gain in the first interval is carried along the trajectory and is reduced
in proportion to every loss. `q` and `dq` are `(parcels, intervals)` as in `adjust_losses`. Returns an array of the
type of `dq`.
"""
function follow_losses(q, dq)
    n, m = size(dq)
    out = zeros(eltype(dq), n, m)
    out[:, 1] = clean.(dq[:, 1])
    @inbounds for t in 2:m, j in 1:n
        d = dq[j, t]
        out[j, t] = clean(d < 0 ? out[j, t-1] * ((q[j, t] + d) / q[j, t]) : out[j, t-1])
    end
    return out
end

"""
    backward_contributions(q, dq, theta, dtheta)

Sources of the parcels in the sink: `(E2P, rows)`, `(E2Q, rows)` and `(H2T, rows)` (`nothing` without `theta`), each a
field `(parcels, intervals)` with the rows of the selected parcels it belongs to.
"""
function backward_contributions(q, dq, theta, dtheta)
    m = size(dq, 2)
    e2q = adjust_losses(view(q, :, 1:m), dq)
    arrival = dq[:, m] .< 0
    e2p = arrival_precipitation(e2q[arrival, :], dq[arrival, m], q[arrival, m])
    h2t = isnothing(theta) ? nothing : (adjust_losses(view(theta, :, 1:m), dtheta), :)
    return (e2p, arrival), (e2q, :), h2t
end

"""
    forward_contributions(q, dq, theta, dtheta)

Fate of the moisture and heat that the parcels gain in the source in the first interval: `(PFE, rows)`, `(QFE, rows)`
and `(TFH, rows)` (`nothing` without `theta`). The precipitation `PFE` is scaled so that it adds up to the first gain.
"""
function forward_contributions(q, dq, theta, dtheta)
    m = size(dq, 2)
    start = dq[:, 1] .> 0
    qfe = follow_losses(q[start, 1:m], dq[start, :])
    pfe = zeros(eltype(qfe), size(qfe))
    pfe[:, 2:end] = qfe[:, 1:end-1] .- qfe[:, 2:end]
    fraction = dq[start, 1] ./ vec(sum(pfe, dims=2))
    pfe = clean.(pfe .* fraction)
    tfh = nothing
    if !isnothing(theta)
        heat_start = dtheta[:, 1] .> 0
        tfh = (follow_losses(theta[heat_start, 1:m], dtheta[heat_start, :]), heat_start)
    end
    return (pfe, start), (qfe, start), tfh
end

"""
    in_region(region, lon, lat)

Which parcels (rows of `lon` and `lat`, two time steps `(parcels, 2)`) are in the region between the two time steps:
their midpoint, rounded to 0.25 degrees, is in the set `region` of `(lon, lat)` tuples (see `read_mask`).
"""
function in_region(region, lon, lat)
    lonmid, latmid = midpoints(lon[:, 2], lon[:, 1], lat[:, 2], lat[:, 1])
    snap(x) = Float64(round(x * 4.0f0) / 4.0f0) + 0.0   # + 0.0 turns -0.0 into 0.0 so that it matches 0.0 in the region
    return [(snap(x), snap(y)) in region for (x, y) in zip(lonmid, latmid)]
end

"""
    attribute(config, region, datesmid, lon, lat, qv, Delta_q, theta, Delta_theta)

Backward (`forward = false`): moisture and heat sources of the parcels that are in the sink region at the end of the
window. Forward (`forward = true`): where the moisture and heat that the parcels gain in the source region at the start
of the window go. The arrays are `(parcels, times)`; `datesmid` are the dates in the middle of the intervals. Returns the
named tuple `(E2P, E2Q, H2T)` (backward) or `(PFE, QFE, TFH)` (forward) of daily `(nx, ny, days)` fields in mm (W/m^2
for heat); the second is `nothing` without `save_q`, the third without `track_heat`.
"""
function attribute(config, region, datesmid, lon, lat, qv, Delta_q, theta, Delta_theta)
    g, flags = grid_definition(config), config["flags"]
    forward, heat = flags["forward"], flags["track_heat"]
    ntimes = size(lon, 2)
    ends = forward ? (1:2) : (ntimes-1:ntimes)
    sel = in_region(region, view(lon, :, ends), view(lat, :, ends))
    contributions = forward ? forward_contributions : backward_contributions
    main, q, h = contributions(qv[sel, :], Delta_q[sel, 2:end], heat ? theta[sel, :] : nothing, heat ? Delta_theta[sel, 2:end] : nothing)

    ls, la = lon[sel, :], lat[sel, :]
    x, y = cells(g, view(ls, :, 2:ntimes), view(ls, :, 1:ntimes-1), view(la, :, 2:ntimes), view(la, :, 1:ntimes-1))

    days = Date.(datesmid)
    daily(a) = stack(d -> dropdims(sum(view(a, :, :, days .== d), dims=3), dims=3), unique(days))
    grid((field, rows)) = daily(ongrid(g.nx, g.ny, view(x, rows, :), view(y, rows, :), field, g.mass))
    names = forward ? (:PFE, :QFE, :TFH) : (:E2P, :E2Q, :H2T)
    return NamedTuple{names}((to_mm(g, grid(main)), flags["save_q"] ? to_mm(g, grid(q)) : nothing, heat ? to_wm2(g, grid(h), 24 * 3600) : nothing))
end

"Parcel ranges of `n` blocks of equal size; the last block takes the remainder."
function parcel_blocks(num_parc, n)
    1 <= n <= num_parc || throw(ArgumentError("the number of blocks must be between 1 and $num_parc, got $n"))
    len = div(num_parc, n)
    return [(i-1)*len+1:(i == n ? num_parc : i * len) for i in 1:n]
end

"Drop the first time step of `window` (a named tuple of `(parcels, times)` arrays) and append the new one, in place."
shift_window(window, new) = map(window, new) do a, b
    copyto!(a, 1, a, size(a, 1) + 1, length(a) - size(a, 1))
    a[:, end] = b
    return a
end

"Sum of the fields of the block results `rs`; `nothing` where a field is off."
combine(rs) = map((fields...) -> isnothing(first(fields)) ? nothing : sum(fields), rs...)

"""
    track(config; blocks=1)

Attribution for every reference date. Backward: the arrival dates from `start_date + dt` to `end_date`, tracked over the
`num_days` before. Forward: the starting dates from `start_date` to `end_date - dt`, tracked over the `num_days` after.
`attribution_{e2p,e2q,h2t}_<date>.nc` (forward `pfe`, `qfe`, `tfh`, dated at the end of the first interval) is written to
`path_data2/attribution`. The first date reads all files; every later date reads only its newest file and shifts the
window.

With MPI initialised (`MPI.Init()`, started with `mpiexec`), every rank takes an equal share of the parcels, the fields are
summed on rank 0 and rank 0 writes them. Within a rank, the parcels are split into `blocks` parts that are processed one
after the other and summed.
"""
function track(config; blocks=1)
    c, d, flags, paths = config["constants"], config["dates"], config["flags"], config["paths"]
    dt, half, forward = c["dt"], Second(c["dt"] * 1800), flags["forward"]
    comm = MPI.Initialized() ? MPI.COMM_WORLD : nothing
    rank, ranks = isnothing(comm) ? (0, 1) : (MPI.Comm_rank(comm), MPI.Comm_size(comm))
    region = Set(read_mask(joinpath(paths["path_mask"], "mask.nc")))
    names = ["lon", "lat", "sh", "Delta_q"]
    flags["track_heat"] && append!(names, ["theta", "Delta_theta"])
    mine = parcel_blocks(c["num_parc"], ranks)[rank+1]
    parcels = [mine[p] for p in parcel_blocks(length(mine), blocks)]
    first_date, last_date = DateTime(d["start_date"], DATEFORMAT), DateTime(d["end_date"], DATEFORMAT)
    references = forward ? (first_date:Hour(dt):last_date-Hour(dt)) : (first_date+Hour(dt):Hour(dt):last_date)

    window = nothing
    for reference in references
        starting, arrival = forward ? (reference, reference + Day(c["num_days"])) : (reference - Day(c["num_days"]), reference)
        dates = isnothing(window) ? (starting:Hour(dt):arrival) : [arrival]
        datesmid = starting+half:Hour(dt):arrival-half
        files = partoutput_files(joinpath(paths["path_data1"], "corrected"), dates)
        new = [read_window(files, p, names) for p in parcels]
        window = isnothing(window) ? new : shift_window.(window, new)

        r = combine([attribute(config, region, datesmid, w.lon, w.lat, w.sh, w.Delta_q, get(w, :theta, nothing), get(w, :Delta_theta, nothing)) for w in window])
        isnothing(comm) || (r = map(field -> isnothing(field) ? nothing : MPI.Reduce(field, +, comm; root=0), r))
        rank == 0 || continue
        tlevel = DateTime.(unique(Date.(datesmid)))
        stamped = forward ? reference + Hour(dt) : reference
        for (name, field) in pairs(r)
            file = joinpath(paths["path_data2"], "attribution", "attribution_$(lowercase(String(name)))_$(stamp(stamped)).nc")
            isnothing(field) || write_attribution(file, config, stamped - half, tlevel, String(name), field)
        end
    end
end
