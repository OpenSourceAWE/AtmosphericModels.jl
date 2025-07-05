```@meta
CurrentModule = AtmosphericModels
```

## Introduction

## Types

### Exported types
```julia
AtmosphericModel
@enum ProfileLaw EXP=1 LOG=2 EXPLOG=3 FAST_EXP=4 FAST_LOG=5 FAST_EXPLOG=6
```

## Functions

### Exported functions
```julia
clear(s::AM)
calc_wind_factor(am::AM, height, profile_law::Int64=am.set.profile_law)
```

```@docs
calc_rho
get_wind
```