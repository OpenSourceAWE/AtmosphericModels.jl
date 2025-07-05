```@meta
CurrentModule = AtmosphericModels
```

## Introduction

## Types

### Exported types
```@docs
ProfileLaw
AtmosphericModel
AtmosphericModel(set::Settings; nowindfield::Bool=false)
```

## Functions

### Wind shear and air density calculation
```@docs
clear
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