# Changelog

## AtmosphericModels v0.3.8 (unreleased)
### Added
- `windfield_path` and `set_windfield_path!`. The `.npz` wind fields are now written to a
  `Scratch.jl` scratchspace instead of `get_data_path()`: they are derived artifacts of ~1.2 GB
  apiece, and writing them next to version-controlled settings made every downstream repo
  accumulate its own copies (7.0 GB in `KiteControllers.jl/data`, 4.7 GB in `V3Kite/data`, ... none
  of them shared) and need its own `.gitignore` rule. The scratchspace is shared by all consumers
  and is removed with the package. `set_windfield_path!(path)` overrides it, `""` restores the
  default. Files in `get_data_path()` are still found and used, so nothing has to be regenerated —
  see `find_windfield`.

### Changed
- `WindField(am, speed)` no longer catches every exception, logs it and returns `nothing`. That
  turned any failure into a `nothing` wind field, which surfaced much later as the
  `wf !== nothing` assertion in `get_wind` with a stack trace nowhere near the cause. It now
  validates the settings up front (`check_windfield_settings`, which names the offending
  `environment.*` key) and lets anything else propagate. Callers testing the result for `nothing`
  need to catch instead.
- The wind field file name now carries `param_digest(set)`, eight hex digits of a SHA-256 over
  `grid_step`, `height_step`, `i_ref`, `alpha`, `avg_height` and `h_ref`. Before, only `grid` and
  the ground wind speed were in the name, so changing any of the other six silently loaded a stale
  file and produced plausible, wrong numbers — the "known limitation" in `docs/src/wind_field.md`,
  whose workaround was deleting every `.npz` by hand. Files under the older names are still found
  and used (`find_windfield`), with a log line saying they cannot be checked against the settings.
  `calc_basename(set)` gained the digest, `grid_basename(set)` is the name without it, and `load`
  lost its `basename` keyword.
- The `.npz` files no longer store the `x`, `y`, `z` coordinate meshgrids, which were half of
  every file (622 MB of 1.24 GB for the default grid) and were used only to derive six scalars
  nothing reads. `grid_axes(am)` rebuilds the axes from `set.grid`/`grid_step`/`height_step`
  instead, so a new file is half the size and loads in half the time. Older files still load, and
  faster than before, because `load` now names the variables it wants instead of reading the whole
  archive. `save` and `load` no longer take or return `x`, `y`, `z`, and `WindField`'s `x`/`y`/`z`
  are the 1D axes rather than 3D meshgrids.
- `use_turbulence` is now applied when the wind field is read (`get_wind`) instead of being baked
  into the stored field, and it is no longer part of the `.npz` filename. One file per ground wind
  speed therefore serves all turbulence intensities, and changing `use_turbulence` no longer
  requires regenerating a ~1.2 GB file. The wind vectors returned are unchanged:
  `sigma = use_turbulence * rel_turbs[idx] * sigma_IEC(v)` as before.
- `WindField` gained a `v_wind_gnd` field, the `set.v_wind_gnds` entry the loaded field was
  generated for, and `load_windfield` returns it as an eighth element. `get_wind` takes the
  `rel_turbs` correction for that speed instead of for `set.v_wind`: the two agreed only because
  `AtmosphericModel(set)` happens to load the field at `set.v_wind`, so any caller constructing
  `WindField(am, speed)` with another speed silently paired one scenario's field with another
  scenario's turbulence intensity.
- **Migration**: existing `windfield_<grid>_1.0_<speed>.npz` files are still valid, just renamed —
  drop the `_1.0` from the name (`load` falls back to the old name and tells you). Files generated
  with any other `use_turbulence` are pre-scaled and should be deleted.

## AtmosphericModels v0.3.7 2026-08-07
### Added
- add `CUSTOM_LOG`, `CUSTOM_EXP` and `CUSTOM_JET` profile laws (`profile_law` 4/5/6), fitting a
  wind profile to `set.heights`/`set.speeds` instead of a fixed `alpha`/`z0`
- add `custom_log`/`custom_exp`, ordinary least squares fits of a logarithmic/power-law profile
- add `custom_jet`, a nonlinear least squares (Levenberg-Marquardt) fit of a power-law
  background plus a superimposed Gaussian jet, `u(z) = c*z^a + U_J*exp(-(z-z_c)^2/(2*sigma^2))`
- cache the `CUSTOM_JET` fit in `AtmosphericModel.jet_cache`, reused while
  `set.heights`/`set.speeds` are unchanged (~200x faster on a cache hit than refitting)
- add examples `plot_custom_exp_log.jl` and `plot_custom_jet.jl`
- add `bench_profile_law.jl`, benchmarking all profile laws, including cold vs. cached `CUSTOM_JET`

### Changed
- Dropped support for Julia 1.10

### Fixed
- fix `new_windfield` generating a different wind field per Julia version: it seeded the global
  RNG (`Random.seed!(1234)`), but `randn`'s array-filling algorithm isn't guaranteed stable across
  Julia versions, which broke the `calc_turbulent_wind` reference values on Julia 1.10. Now uses a
  `StableRNG(1234)` instead, so the wind field is byte-identical everywhere.

## AtmosphericModels v0.3.6 2026-08-05
### Added
- add `calc_turbulent_wind`, moved here from `KiteModels` (kite/tether wind vectors in the ENU
  frame, built on top of `get_wind`)
- add `bin/release`, which posts the latest `CHANGELOG.md` release notes to
  `OpenSourceAWE/AtmosphericModels.jl` issue #1 to trigger `JuliaRegistrator`, after checking that
  the working tree is clean, `Project.toml`'s version matches the changelog entry, and
  `Manifest-v1.12.toml.default` matches `Manifest-v1.12.toml`

### Fixed
- fix `get_wind` docstring: the `upwind_dir` default is `-π/4`, not `0.0`
- fix `get_wind` to detect which of the field's first two dimensions is the long (along-wind) one
  at runtime, instead of assuming dimension 1 is always long. That assumption broke any
  `set.grid` with the short dimension first (e.g. KiteUtils' own default `[100, 4050, ...]`,
  as opposed to this package's own `data/settings.yaml` default `[4050, 100, ...]`).

## AtmosphericModels v0.3.5 2026-05-30
### Added
- add `CONSTANT` profile law (no wind shear, `profile_law = 0`)
- add example `plot_windshear_zero.jl`

### Changed
- allow ControlPlots 0.3 in examples compat

### Fixed
- fix typo in `bin/run_julia`: `JULIA_PKG_SERVER_REGISTRY_PREFERANCE` → `JULIA_PKG_SERVER_REGISTRY_PREFERENCE`

## AtmosphericModels v0.3.4 2026-05-03
### Added
- add `upwind_dir` to `get_wind`
- add helper scripts `bin/install`, `bin/setup_env`, and `bin/jetls`
- add `menu()` to `bin/run_julia`
- add `.markdownlint.json`

### Changed
- improve `bin/run_julia`
- update default manifests for Julia 1.11 and 1.12
- update documentation and README
- remove TestEnv usage from project tooling

### Fixed
- fix installation script behavior
- fix warning in windfield code

## AtmosphericModels v0.3.3 2026-03-17
### Added
- the files `.zenodo.json` and `CITATION.cff`
### Changed
- support Julia 1.12

## AtmosphericModels v0.3.2 2025-08-26
### Added
- add KiteUtils 0.11 compat

## AtmosphericModels v0.3.1 2025-07-14
### Fixed
- when calculating the filename for the windfield to load, the `rel_sigma` parameter was ignored
### Changed
- better error message if loading the windfield fails

## AtmosphericModels v0.3.0 2025-07-08
### Changed
- BREAKING: When constructing an atmospheric model, you MUST pass the parameter set::Settings. This ensures that all parts of the simulation use the same settings struct, and that you can run different simulations with different settings in parallel.
- removed FAST_EXP, FAST_LOG and FAST_EXPLOG because they were error prone (did not deliver the correct result when changing settings.yaml)

### Added
- The function `get_wind(am, x, y, z, t)` which returns a wind vector for the given position and time. It creates a 3D wind field if it does not exist in the data folder. The parameters of this wind field are configured in `settings.yaml`.
- Documenter generated documentation.
- all files have now a license attached. You can check that with `pipx run reuse lint`.
- many examples
- a GUI to investigate the 3D wind field
