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