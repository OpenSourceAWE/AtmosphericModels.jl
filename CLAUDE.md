<!--
SPDX-FileCopyrightText: 2026 Uwe Fechner
SPDX-License-Identifier: MIT
-->

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this package is

AtmosphericModels.jl provides the wind/atmosphere model used by the "Julia Kite Power Tools" family
(sibling packages: `KiteModels`, `KiteUtils`, `WinchModels`, `KitePodModels`, all under the
OpenSourceAWE / aenarete GitHub orgs). It has two responsibilities:

1. **Deterministic profiles**: air density vs. height (`calc_rho`) and mean wind speed vs. height
   under a chosen wind profile law (`calc_wind_factor`).
2. **Turbulent 3D wind fields**: generates/loads a Mann-model turbulence field (`WindField`) and
   returns a wind vector at any `(x, y, z, t)` via `get_wind`, using Taylor's frozen-turbulence
   hypothesis (the field is advected along the mean wind direction rather than recomputed per step).
   `KiteModels.jl` calls into this for turbulent-wind simulations.

The Mann-model turbulence code in `src/windfield.jl` is a Julia port of a Matlab module (`Mann.m` by
René Bos), itself based on Mann (1994) and Mann (1998) — see the file's module docstring for the
paper citations. When touching `create_windfield`/`create_grid`, keep in mind the code intentionally
mirrors Python/Matlab array conventions in a few places (e.g. `create_grid` returns `(Y, X, Z)`, not
`(X, Y, Z)`, "to match the Python (y, x, z) order" — the comment in the code is load-bearing, not a
mistake to "fix").

## Architecture

### Single small module, two files

`src/AtmosphericModels.jl` defines the module, the `AtmosphericModel`/`WindField` structs, the
`ProfileLaw` enum (`CONSTANT`/`EXP`/`LOG`/`EXPLOG`), `calc_rho`, and `calc_wind_factor`, then
`include`s `src/windfield.jl` for everything turbulence-related (`get_wind`, `new_windfield`,
`new_windfields`, `create_windfield`, `create_grid`, `load`/`save` of `.npz` wind-field files).
There is no submodule split — both files contribute to the same `AtmosphericModels` module.

- `AtmosphericModel(set::Settings; nowindfield=false)` — the constructor. If
  `set.use_turbulence > 0` and `nowindfield=false`, it eagerly loads (or generates, if missing) a
  `WindField` for `set.v_wind`. Precomputes `rho_zero_temp` from `set.temp_ref`/`set.rho_0`; call
  `clear(am)` after mutating `set.temp_ref` to recompute it (see `test/runtests.jl` for the pattern).
- `WindField` — an immutable `@kwdef struct` holding the turbulence grid (`u`,`v`,`w` 3D arrays,
  `x`,`y`,`z` axes) plus a custom `Base.getproperty` that synthesizes `x_range`/`y_range` from the
  stored min/max fields.
- `calc_wind_factor(am, height, profile_law)` dispatches on `Val{Int(profile_law)}` for the hot path
  (`calc_wind_factor1/2/3` per law); there's also a non-`Val` `Int64`-dispatched overload with a
  runtime `if`/`elseif` for callers that don't know the law at compile time.
- Wind field `.npz` files live in `data/` and are named via `calc_basename(set)` (from
  `set.grid`) + ground wind speed + `rel_sigma` — see `calc_full_name`/`load_windfield`. If a needed
  file is missing, `load` transparently calls `new_windfield` to generate it (can take ~30s+;
  `Random.seed!(1234)` makes generation deterministic).

### Configuration

All physical parameters live under the `environment:` section of `data/settings.yaml` (see
`docs/src/settings.md` for annotated examples for two real sites, Cabauw and Maasvlakte NL).
`system.yaml` selects which settings file is active. Turbulence-specific keys (`v_wind_gnds`,
`avg_height`, `rel_turbs`, `i_ref`, `v_ref`, `grid`, `height_step`, `grid_step`) are only needed when
`use_turbulence > 0`; with `use_turbulence == 0` no wind field is loaded at all. `Settings` objects
come from `KiteUtils`, not this package — load with `set_data_path(...)` +
`load_settings("system.yaml"; relax=true)`.

## Development commands

Workspace-based like the sibling packages: `Project.toml` declares
`[workspace] projects = ["examples", "test", "docs"]`.

- **Install/setup**: `cd bin && ./install` (juliaup + `setup_env`); `./install --update` to refresh.
- **Launch a dev REPL**: `./bin/run_julia` (defaults to `using KiteUtils: menu`; forwards script args
  if given).
- **Run the test suite** (project = `test/`):
  Always use the ex tool of Kaimon to run the tests:
  ```julia
  include("test/runtests.jl")
  ```
  Unlike `KiteModels.jl`, this is a flat script (no directory walk over `test-*.jl` files) — it
  loads `data/system.yaml`, builds one shared `am`, includes `test_windfield.jl`, then runs the
  `calc_wind_factor`/`calc_rho` testsets inline. To run just the windfield tests, include
  `test/test_windfield.jl` directly (after setting up `am` as `runtests.jl` does).
- **Run examples**: `include("examples/menu.jl")` for the interactive menu (options include
  `get_wind`, `plot_windshear*`, `new_windfields` to regenerate all turbulence files for
  `set.v_wind_gnds`, and `delete_windfields` to remove cached `.npz` files from `data/`).
- **Build docs locally**:
  ```julia
  using Pkg; Pkg.activate("docs"); include("docs/make.jl"); Pkg.activate(".")
  ```
- **License linting**: every file must carry an SPDX header (see `REUSE.toml` for path-based
  annotation defaults); run via `pipx run reuse lint`.

## Coding style

Same conventions as the other Julia Kite Power Tools packages (see `KiteModels.jl`'s `CLAUDE.md` /
`docs/src/advanced.md` for the full list): 120-char line limit, named constants over magic numbers,
`\cdot` for dot products, space around binary operators, install `Revise` globally (never as a
project dependency).
