```@meta
CurrentModule = AtmosphericModels
```

## Introduction
Most functions need an instance of the struct `AtmosphericModel` as first parameter, which can be created using the following code:
```julia
using AtmosphericModels, KiteUtils

set_data_path("data")
set = load_settings("system.yaml"; relax=true)
am::AtmosphericModel = AtmosphericModel(set)
```
This requires that the files `system.yaml` and `settings.yaml` exist in the folder `data`. See also [Settings](@ref). The parameter `relax=true` allows loading a yaml file that does not contain all sections needed to run a kite power system simulation. This is useful if you want to use this package for
other purposes than simulating kite power systems.

## Types

### Exported types
```@docs
ProfileLaw
AtmosphericModel
AtmosphericModel(set::Settings; nowindfield::Bool=false)
```
### Private types
```@docs
WindField
```

## Functions

### Wind shear and air density calculation
```@docs
calc_rho
calc_wind_factor
```
### Wind turbulence calculation
```@docs
get_wind
rel_turbo
new_windfield
new_windfields
```