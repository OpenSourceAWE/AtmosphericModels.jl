## Move function calc_turbulent_wind from KiteModels to this package

"""
    calc_turbulent_wind(am, pos, upwind_dir, t)

Calculate the wind velocity vectors at the kite and at the mid-tether point.

When `am.set.use_turbulence == 0.0`, returns the mean wind based on a log/power-law
height profile. When `use_turbulence != 0.0`, returns the fully turbulent wind vectors
looked up from the pre-computed wind field via `get_wind`.

Parameters:
- `am`:         atmospheric model (settings are read from `am.set`)
- `pos`:        3D position of the kite [m]; `pos[3]` is used as height (clamped to 6 m minimum)
- `upwind_dir`: upwind direction in radians; zero is north, clockwise positive
- `t`:          current simulation time [s]

Returns a tuple `(v_wind, v_wind_tether)` where:
- `v_wind`:        wind velocity vector at kite height [m/s]
- `v_wind_tether`: wind velocity vector at half the kite height [m/s]
"""
function calc_turbulent_wind(am, pos, upwind_dir, t)