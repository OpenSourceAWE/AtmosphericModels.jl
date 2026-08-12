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
mistake to "fix"). Both build on `grid_axes(am)`, the single definition of the grid's coordinate
axes; the `.npz` stores only `u`, `v`, `w` and `param`, so the axes are rebuilt, never read back.

## Architecture

### Single small module, two files

`src/AtmosphericModels.jl` defines the module, the `AtmosphericModel`/`WindField` structs, the
`ProfileLaw` enum (`CONSTANT=0`/`EXP=1`/`LOG=2`/`EXPLOG=3`/`CUSTOM_LOG=4`/`CUSTOM_EXP=5`/
`CUSTOM_JET=6`), `calc_rho`, `calc_wind_factor`, and the custom-profile fits
(`custom_log`/`custom_exp`/`custom_jet`), then `include`s `src/windfield.jl` for everything
turbulence-related (`get_wind`, `get_wind!`, `calc_turbulent_wind`, `new_windfield`,
`new_windfields`, `create_windfield`, `create_grid`, `load`/`save` of `.npz` wind-field files).
There is no submodule split — both files contribute to the same `AtmosphericModels` module.
`calc_turbulent_wind(am, pos, t; upwind_dir)` — which returns the wind vector at the kite plus the
one at half its height for the tether — moved here from `KiteModels.jl` in Feb 2026; `KiteModels`'
`set_v_wind_ground!` is its only caller in the family.

- `AtmosphericModel(set::Settings; nowindfield=false)` — the constructor. If
  `set.use_turbulence > 0` and `nowindfield=false`, it eagerly loads (or generates, if missing) a
  `WindField` for `set.v_wind`. That throws if it cannot: `check_windfield_settings` validates the
  `environment.*` keys a field needs first, and nothing is swallowed. Precomputes `rho_zero_temp` from `set.temp_ref`/`set.rho_0`; call
  `clear(am)` after mutating `set.temp_ref` to recompute it (see `test/runtests.jl` for the pattern).
- `WindField` — an immutable `@kwdef struct` holding the turbulence grid (`u`,`v`,`w` 3D arrays,
  `x`,`y`,`z` axes) plus a custom `Base.getproperty` that synthesizes `x_range`/`y_range` from the
  stored min/max fields. Its `v_wind_gnd` records which `set.v_wind_gnds` entry the field was
  generated for; `get_wind` reads it back to pick the matching `rel_turbs`, so a field loaded for
  a speed other than `set.v_wind` keeps its own turbulence intensity.
- `calc_wind_factor(am, height, profile_law)` dispatches on `Val{Int(profile_law)}` for the hot path
  (`calc_wind_factor1` … `calc_wind_factor6`, one per non-`CONSTANT` law); there's also a non-`Val`
  `Int64`-dispatched overload with a runtime `if`/`elseif` for callers that don't know the law at
  compile time.
- The three `CUSTOM_*` laws fit a profile to the `heights`/`speeds` pairs in the settings instead of
  using `alpha`/`z0`: `custom_log` (least squares on `u = a·log z + b`), `custom_exp` (least squares
  in log-log space, `u = c·z^a`), and `custom_jet` (a 5-parameter power-law-plus-Gaussian
  low-level-jet model, fitted by the hand-rolled Levenberg-Marquardt in `fit_jet`). All divide by
  `set.v_wind` to return a factor. `fit_jet` is the expensive one, so `jet_coeffs(am)` memoizes it in
  `am.jet_cache`, which stores copies of `heights`/`speeds` and compares them by value — so an
  in-place edit of the settings does invalidate the cache, and `clear(am)` (which only recomputes
  `rho_zero_temp`) does not need to touch it.
- Wind field `.npz` files live in a `Scratch.jl` scratchspace (`windfield_path()`, redirectable
  with `set_windfield_path!`) and are named via `calc_basename(set)` (`set.grid` + a digest of
  every other setting the field depends on) + ground wind speed — see
  `calc_full_name`/`param_digest`/`find_windfield`/`load_windfield`. Files
  written before v0.3.8 sat in `get_data_path()` under a name without the digest; `find_windfield`
  still finds them, preferring the digested name. They store the field at
  the reference intensity; `set.use_turbulence` and `rel_turbs` are applied in `get_wind`, so one
  file serves every intensity. If a needed file is missing, `load` transparently calls
  `new_windfield` to generate it (can take ~30s+; the `StableRNG(1234)` seed makes generation
  deterministic), falling back to the pre-0.3.8 `..._1.0_<speed>.npz` name if that file exists.

### Configuration

All physical parameters live under the `environment:` section of `data/settings.yaml` (see
`docs/src/settings.md` for annotated examples for two real sites, Cabauw and Maasvlakte NL — the
latter is `data/settings_nearshore.yaml`, selected by `data/system_nearshore.yaml`).
`system.yaml` selects which settings file is active. `heights`/`speeds` are only read by the
`CUSTOM_*` profile laws. Turbulence-specific keys (`v_wind_gnds`,
`avg_height`, `rel_turbs`, `i_ref`, `v_ref`, `grid`, `height_step`, `grid_step`) are only needed when
`use_turbulence > 0`; with `use_turbulence == 0` no wind field is loaded at all. `Settings` objects
come from `KiteUtils`, not this package — load with `set_data_path(...)` +
`load_settings("system.yaml"; relax=true)`.

## Development commands

Workspace-based like the sibling packages: `Project.toml` declares
`[workspace] projects = ["examples", "test", "docs"]`.

- **Install/setup**: `cd bin && ./install` (juliaup). It takes no options — unlike `KiteModels.jl`'s
  installer there is no `--update` flag; just re-run it.
- **Launch a dev REPL**: `./bin/run_julia` (defaults to `using KiteUtils: menu`; forwards script args
  if given).
- **Run the test suite** (project = `test/`):
  Always use the ex tool of Kaimon to run the tests:
  ```julia
  include("test/runtests.jl")
  ```
  Unlike `KiteModels.jl`, this is a flat script (no directory walk over `test-*.jl` files) — it
  loads `data/system.yaml`, builds one shared `am`, includes `test_windfield.jl` and
  `test_custom_profiles.jl`, then runs the `calc_wind_factor`/`calc_rho` testsets inline. To run just
  the windfield tests, include `test/test_windfield.jl` directly (after setting up `am` as
  `runtests.jl` does).
- **Run examples**: `include("examples/menu.jl")` for the interactive menu (options include
  `get_wind`, `plot_windshear*`, `plot_custom_exp_log`/`plot_custom_jet` for the `CUSTOM_*` laws,
  `new_windfields` to regenerate all turbulence files for `set.v_wind_gnds`, and
  `delete_windfields`, which deletes every `.npz` in `windfield_path()` — the scratchspace, not
  `data/`).
- **Build docs locally**:
  ```julia
  using Pkg; Pkg.activate("docs"); include("docs/make.jl"); Pkg.activate(".")
  ```
- **License linting**: every file must carry an SPDX header (see `REUSE.toml` for path-based
  annotation defaults); run via `bin/reuse_lint`, which tries `pipx run reuse lint` and falls back to
  a plain `reuse`/`python3 -m reuse` if pipx is missing or broken.
- **Cut a release**: `./bin/release` — takes no arguments, must be run on `main`, extracts the top
  section of `CHANGELOG.md` as release notes and posts them to GitHub issue #1 (the JuliaRegistrator
  trigger issue) via `gh`.

## Coding style

Same conventions as the other Julia Kite Power Tools packages (see `KiteModels.jl`'s `CLAUDE.md` /
`docs/src/advanced.md` for the full list): 120-char line limit, named constants over magic numbers,
`\cdot` for dot products, space around binary operators, install `Revise` globally (never as a
project dependency).
