# AtmosphericModels wind-field issues

Findings from debugging a `use_turbulence: 0.2` run of `examples/simple_fig8.jl` on
2026-08-08. All line numbers refer to `AtmosphericModels/src/windfield.jl` in the local
checkout at `~/repos/AtmosphericModels` (v0.3.7); v0.3.6 as installed is identical in
every line quoted here. Nothing below is a bug in SimpleKiteControllers.jl — this file
is a bug report to file upstream.

## 1. `use_turbulence` is baked into the stored field instead of applied at lookup

`new_windfield` scales the Mann field before saving:

```julia
sigma1 = am.set.use_turbulence * calc_sigma1(am, v_wind_gnd)   # :447
```

and `calc_full_name` puts `use_turbulence` into the filename (`:85`), so every intensity
gets its own 1.2 GB file. But `use_turbulence` is a **pure linear scale factor** on a
zero-mean field: the field for 0.5 is exactly `0.5 *` the field for 1.0, same RNG seed,
same everything. Applying it at lookup time in `get_wind` — where `rel_turbs[idx]` is
already applied the same way (`:381-383`) — would make one file per wind speed serve all
intensities, and would make changing the intensity free instead of a 1.2 GB regeneration.

Measured cost of not doing this, on this machine today: 14 files, 18 GB total, of which
five are the same 8.2 m/s field at `use_turbulence` ∈ {0.1, 0.4, 0.5, 0.7, 1.0} and four
are the same 9.5 m/s field at {0.1, 0.5, 0.9, 1.0}. Seven of the fourteen files are
redundant — ~8.7 GB — and would collapse to zero under a lookup-time scale.

Note this is *not* a claim that multiplying by both `use_turbulence` and `rel_turbs` is
wrong. Per `docs/src/wind_field.md`, `rel_turbs` corrects the IEC sigma down to the
measured Cabauw values and `use_turbulence` scales relative to Cabauw, so
`sigma = use_turbulence * rel_turbs[idx] * sigma_IEC(v)` is the intended composition and
`use_turbulence: 1.0` correctly reproduces Cabauw. The problem is only *where* the first
factor is applied.

## 2. `rel_turbo` is indexed by `set.v_wind`, not by the speed the field was built for

The field is chosen by the argument to `WindField(am, speed)`:

```julia
idx = findmin(abs.(am.set.v_wind_gnds .- speed))[2]   # load_windfield, :119
```

but the correction factor is chosen independently, from the settings:

```julia
rel_turbo(am::AtmosphericModel, v_wind = am.set.v_wind) # :66 — get_wind calls it as rel_turbo(am), :325
```

The two agree only because `AtmosphericModel(set)` happens to call
`WindField(am, am.set.v_wind)`. Any other caller — and V3Kite does construct wind fields
at a speed of its own choosing during settling — silently pairs the field generated for
one scenario with the correction factor of another. `rel_turb` should come from the index
the loaded field was actually built with; storing that index (or the speed) on the
`WindField` and reading it back in `get_wind` is the minimal fix.

## 3. The `.npz` files belong in a scratchspace, not in `get_data_path()`

```julia
path = get_data_path() * "/"   # calc_full_name, :86
```

These are derived artifacts: deterministic given the settings plus a hardcoded seed, and
regenerable in ~30-60 s. They are written straight into whatever data directory happens
to be active, i.e. next to version-controlled YAML settings, at 1.2 GB apiece. Every
downstream repo therefore accumulates its own copies and needs its own `.gitignore` rule
— today `KiteControllers.jl/data` holds 7.0 GB, `V3Kite/data` 4.7 GB,
`KiteModels.jl/data` 3.5 GB, `AtmosphericModels/data` 1.2 GB, none of them shared.

`Scratch.jl` (`@get_scratch!("windfields")`) is the right home, exactly as V3Kite already
does for its model and settled-geometry binaries. That gets a single shared cache across
all consumers, automatic cleanup when the package is removed, and no risk of a 1.2 GB
file being committed. It also makes issue 5 below fixable, since a scratchspace path can
carry a hash of every relevant parameter without the filename becoming unreadable.

## 4. Half of every file is a coordinate meshgrid that is never used

`create_grid` returns full 3D arrays from `ndgrid` (`:144-158`) and `save` writes `x`,
`y`, `z` alongside `u`, `v`, `w` (`:96-104`). For the default grid that is
51 × 2026 × 251 = 25.9 M points, so three redundant `Float64` arrays of 207 MB each:
**622 MB of the 1,244,872,856-byte file is reconstructible from `grid`, `grid_step` and
`height_step`.**

They are not used for anything either. On load they serve only to compute six scalars
(`:31-36`), and those six fields (`x_max`, `x_min`, …) are read nowhere in the package —
`get_wind` derives its indices arithmetically from `grid_step`/`height_step`, not from
the stored coordinates. Dropping `x`, `y`, `z` from the file halves it and halves load
time and peak RAM.

## 5. The filename covers only 3 of the ~9 parameters that determine the field

Already flagged as a "Known limitation" at the end of `docs/src/wind_field.md`, with
"delete all `*.npz` files by hand" as the workaround. `height_step`, `grid_step`, `i_ref`,
`alpha`, `avg_height`, `h_ref`, `profile_law` and `z0` all change the field but not its
name, so a stale file is loaded silently and the run produces plausible, wrong numbers.
A content hash over the full parameter set turns this from a documented footgun into a
non-issue, and is natural to do at the same time as issue 3.

## 6. A failure to build the field is swallowed and resurfaces much later

```julia
catch e
    @error "Error reading wind field!"
    println("Caught exception: ", e)
    return nothing          # :38-42
end
```

`WindField` catches everything, logs, and returns `nothing`. Construction then completes
normally, and the run dies later, elsewhere, on

```julia
@assert wf !== nothing "Wind field is not initialized"   # get_wind, :320
```

with a stack trace pointing at `apply_turbulence!` — nowhere near the real cause. This is
what made today's failure expensive to diagnose: the actual exception was
`BoundsError(Int64[], (1,))` from `calc_basename` (`:456`) indexing an empty
`set.grid`, because `data/settings_fig8_*.yaml` did not declare `grid` and
`KiteUtils.Settings` defaults it to `Int64[]`. Rethrowing, or at minimum validating
`set.grid`, `set.grid_step` and `set.height_step` in the `AtmosphericModel` constructor
with a message naming the missing key, would have made it a one-line diagnosis.

(The SimpleKiteControllers side of this is fixed: `grid: [100, 4050, 500, 70]` is now
declared in all three `data/settings_fig8_*.yaml`.)

## 7. `new_windfield` reseeds the global RNG

```julia
Random.seed!(1234)   # :445
```

This perturbs the caller's random stream as a side effect of a cache miss, so a
simulation is not reproducible across a run that generated a field and a run that found
one on disk. It should seed a local `Xoshiro(1234)` and thread it through
`create_windfield` instead.
