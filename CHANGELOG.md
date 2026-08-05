# Changelog

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
