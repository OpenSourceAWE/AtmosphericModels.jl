## Unreleased

## Changed
- BREAKING: When constructing an atmospheric model, you MUST pass the parameter set::Settings. This ensures that all parts of the simulation use the same settings struct, and that you can run different simulations with different settings in parallel.

## Added
- The function `get_wind(am, x, y, z, t)` which returns a wind vector for the given position and time. It creates a 3D wind field if it does not exist in the data folder. The parameters of this wind field are configured in `settings.yaml`.