### Added
- `plot(am::AtmosphericModel; threshold)` in a new Makie package extension, loaded with
  `MakieControlPlots`: see-through 3D isosurfaces of the gusts (red) and lulls (blue) of the
  along-wind turbulence beyond ±`threshold`, shown in the README.
- `plot_interactive(am)`: the same plot with a slider for the threshold.
