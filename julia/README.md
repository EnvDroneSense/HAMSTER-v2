# HAMSTER.jl

Julia port of the Python code in `../python`. The Python code is the reference.

| Path | Holds |
|---|---|
| `Project.toml` | package definition and version bounds |
| `namelist.toml` | settings; replace every `set_...` placeholder before a run |
| `src/HAMSTER.jl` | the module |
| `src/utilities.jl` | midpoints, potential temperature, gridding |
| `src/io.jl` | namelist, sink mask, partoutput and ERA5 reading, output files |
| `src/preprocessing.jl` | moisture and heat changes, bias correction |
| `src/tracking.jl` | WaterSip attribution, backward and forward |
| `scripts/run.jl` | start script |

```
julia --project=julia julia/scripts/run.jl preprocess julia/namelist.toml
julia --project=julia julia/scripts/run.jl track julia/namelist.toml [blocks]
mpiexec -n 32 julia --project=julia julia/scripts/run.jl track julia/namelist.toml [blocks]
```

Tracking splits the parcels over the MPI ranks; rank 0 sums the fields and writes them. `blocks` splits the parcels of
each rank into that many parts, processed one after the other. On a cluster, point MPI.jl at the system MPI first
(`MPIPreferences.use_system_binary()`). Each run logs its total time and peak memory when it finishes.

## Output

| Step | Reads | Writes |
|---|---|---|
| `preprocess` | `path_data1/partoutputs/partoutput_<date>.nc` (raw) | `path_data1/corrected/partoutputs/partoutput_<date>.nc`; with `save_diagnose` also `path_data1/diagnosis/diagnosis_<date>.nc` |
| `track` | `path_data1/corrected/partoutputs/` | `path_data2/attribution/attribution_<e2p,e2q,h2t>_<date>.nc` (forward: `pfe`, `qfe`, `tfh`) |

Unlike the Python code, the raw partoutput files are never overwritten or deleted, so the preprocessing can be re-run.
The corrected files have the layout of the raw files (`lon`, `lat`, `sh`, `Delta_q` and, with `track_heat`, `theta` and
`Delta_theta`), so the corrected directory is the only place that differs from the Python output. Diagnosis and
attribution files have the same dimensions and variables as the Python ones.

## Known bugs kept for parity

These Python behaviours are reproduced on purpose, so that the Julia results can be compared with the Python results.

| # | Priority | Script | Behaviour kept |
|---|---|---|---|
| F6 | medium | `tracking.jl` (`adjust_losses`) | The loss fraction `(q + dq) / q` mixes the corrected `Delta_q` with the raw `q`; it can be negative after an upward correction, which zeroes all earlier contributions |
| F8 | medium | `tracking.jl` (`adjust_losses`, `follow_losses` with `theta`) | The heat loss factor uses the absolute theta (about 300 K), so losses barely discount earlier gains; H2T has no arrival filter or normalisation |
| F9 | medium | `tracking.jl` (`forward_contributions`, `attribute`) | Forward QFE is a moisture state summed over all time steps of a day (8 times a daily mean at `dt` = 3 h) |
| F10 | medium | `tracking.jl` (`arrival_precipitation`) | The `dq / q` factor cancels: 100 % of the arrival precipitation is attributed, and parcels with a zero sum of E2Q are dropped |
| F12 | low | `io.jl` (`read_mask`), `tracking.jl` (`in_region`) | A midpoint that rounds to lon 180.0 is never in the sink, because the mask grid ends at 179.75; no wrap at the dateline |
| F14 | low | `preprocessing.jl` (`delta`) | The 1e-4 threshold for a change to count is absolute, the same for kg/kg and K, and does not scale with `dt` |
| F19 | low | `io.jl` (`read_window`) | NaN latitudes become -1, so those parcels are counted at latitude -1 instead of being skipped |
