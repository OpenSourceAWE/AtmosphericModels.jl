```@meta
CurrentModule = AtmosphericModels
```

## Settings
The parameters af the atmospheric model can be configured in the section `environment` of the `settings.yaml` file in the `data` folder.

The file `system.yaml` specifies which `yaml` files are used to configure
the current project.

### Example for system.yaml
```yaml
system:
    project: "settings.yaml"  # simulator settings
```
Often additional `yaml` files, for example for the controller settings are used.

### Example for settings.yaml
```yaml
environment:
    v_wind: 5.324            # wind speed at reference height          [m/s]
    upwind_dir: -90.0        # upwind direction                        [deg]
    temp_ref: 15.0           # temperature at reference height         [°C]
    height_gnd: 0.0          # height of groundstation above see level [m]
    h_ref:  6.0              # reference height for the wind speed     [m]

    rho_0:  1.225            # air density at zero height and 15 °C    [kg/m³]
    alpha:  0.234            # exponent of the wind profile law
    z0:     0.0002           # surface roughness                       [m]
    profile_law: 1           # 1=EXP, 2=LOG, 3=EXPLOG, 4=FAST_EXP, 5=FAST_LOG, 6=FAST_EXPLOG
    # the following parameters are for calculating the turbulent wind field using the Mann model
    use_turbulence: 1.0      # turbulence intensity relative to Cabauw, NL
    v_wind_gnds: [3.483, 5.324, 8.163] # wind speeds at ref height for calculating the turbulent wind field [m/s]
    avg_height: 200.0        # average height during reel out          [m]
    rel_turbs:   [0.342, 0.465, 0.583] # relative turbulence at the v_wind_gnds
    i_ref: 0.14              # is the expected value of the turbulence intensity at 15 m/s.
    v_ref: 42.9              # five times the average wind speed in m/s at hub height over the full year    [m/s]
                             # Cabauw: 8.5863 m/s * 5.0 = 42.9 m/s
    grid: [4050, 100, 500, 70] # grid size nx, ny, nz and minimal height z_min                              [m]
    height_step: 2.0         # use a grid with 2m resolution in z direction                                 [m]
    grid_step:   2.0         # grid resolution in x and y direction    [m]
```

## Remarks
- If the parameter `use_turbulence` is zero, no windfield is loaded.