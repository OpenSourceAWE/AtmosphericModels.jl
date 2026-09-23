# Changelog

## AtmosphericModels v0.3.12 (unreleased)
### Added
- `AMSettings`, the `environment:` fields of `KiteUtils.Settings` on their own, and
  `AMSettings(file)` to read them from the `environment:` section of any yaml file, applying
  `use_wind_vec` and rejecting an invalid `profile_law` as `load_settings` does.
  `AtmosphericModel` accepts an `AMSettings` as well as a `Settings`, so packages that only need
  the atmosphere no longer need a `system.yaml` and KiteUtils' project layout.

## AtmosphericModels v0.3.11 - 2026-09-24
### Changed
- Dropped support for Julia 1.11. CI tests 1.12 and 1.13, `Manifest-v1.11.toml.default` is no
  longer tracked, `bin/install` offers 1.12 and 1.13, and the README and docs ask for Julia 1.12 or
  later.
- `KiteUtils` compat widens to `"0.12, 0.13"`.
- `bin/install` takes `-y` (run without a terminal), `--update` (update the live manifest and leave
  the tracked `.default` alone) and `-h`. It installs the tracked manifest rather than re-resolving
  it, and no longer changes the juliaup default or appends a `jl` alias to the shell profile.

## AtmosphericModels v0.3.10 - 2026-09-11
### Added
- Support Julia 1.13: `bin/install` offers it as a version choice and accepts it as the detected
  version, `Manifest-v1.13.toml.default` is the tracked default manifest for it, and
  `Manifest-v1.13.toml` (already gitignored) now also ignores the `paper/paper.pdf` build output.
  `Project.toml` compat bounds for `LinearAlgebra`, `Printf`, `Random`, `Statistics` and `julia`
  gained `1.13`, and `SHA` widened to `0.7, 1.0`.
- Add a JOSS paper draft (`paper/paper.md`, `paper/paper.bib`, `paper/wind_profile.png`,
  `paper/build`) describing the package, and license it under CC-BY-4.0 (`LICENSES/CC-BY-4.0.txt`,
  annotated for `paper/*.*` in `REUSE.toml`).
- Add tests for the `Int64`-dispatched `calc_wind_factor` overload: each branch is checked against
  its `Val`-dispatched counterpart at two heights, the default `profile_law` argument is checked to
  fall back to `am.set.profile_law`, and an out-of-range `profile_law` is checked to raise a
  `DomainError`.

### Fixed
- fix `docs/Project.toml`: add a `[sources]` entry pointing `AtmosphericModels` at `..`, so the docs
  environment resolves the in-repo package instead of a registered release.

## AtmosphericModels v0.3.9 - 2026-08-12
### Changed
- Bump `KiteUtils` to 0.12.

## AtmosphericModels v0.3.8 - 2026-08-09
### Added
- The `interpolate` keyword of `get_wind` is implemented. It used to be a documented keyword that
  returned `nothing` (a `TODO` left over from the Python original, which used
  `ndimage.map_coordinates`). With `interpolate=true` the turbulence is now interpolated trilinearly
  between the eight surrounding grid points instead of read at the nearest one, which removes the
  steps a kite flying through the field sees, at about 1.8x the cost of the lookup (29 ns → 51 ns
  per position for the vector method). The horizontal axes wrap with the same period the
  nearest-grid-point lookup uses, the vertical one is clamped at the top layer. The keyword is
  available on all `get_wind` methods and on `calc_turbulent_wind`; the default stays `false`, so
  nothing changes for existing callers.
- `get_wind(am, positions, t; upwind_dir)`, a method that takes a vector of 3D positions and returns
  a `Vector{SVec3}` of wind vectors, plus the in-place `get_wind!(res, am, positions, t; upwind_dir)`.
  Everything that does not depend on the position — `rel_turbo` (which allocates), the sine/cosine of
  the wind direction and the `am.set` lookups — is computed once per call instead of once per
  position, which makes it about 2.4x faster than the scalar `get_wind` in a loop (100 positions:
  7.2 µs and 203 allocations → 3.0 µs and 5, or 2 with `get_wind!`). The results are bit-identical
  to the scalar method; the scalar path is unchanged in speed.
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
