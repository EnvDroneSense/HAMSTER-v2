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
| `test/runtests.jl` | tests |

```
julia --project=julia julia/scripts/run.jl preprocess julia/namelist.toml
julia --project=julia julia/scripts/run.jl track julia/namelist.toml [blocks]
julia --project=julia julia/test/runtests.jl
```

`blocks` splits the parcels into that many parts, processed one after the other. The tests need the reference files in
`test/reference/`; they are not in git and are written by `make_reference.py`.
