# HAMSTER.jl

Julia port of the Python code in `../python`, started from an empty package. The Python code is the reference.

| Path | Holds |
|---|---|
| `Project.toml` | package definition (standard library only so far) |
| `src/HAMSTER.jl` | the module; ported parts are added here |
| `scripts/run.jl` | start script: `preprocess` or `track` with a namelist |
| `test/runtests.jl` | tests |

```
julia --project=julia julia/scripts/run.jl track python/namelist.input
julia --project=julia julia/test/runtests.jl
```

Packages are added with `Pkg.add` when a part needs them.
