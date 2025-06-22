# AtmosphericModels

[![Build Status](https://github.com/OpenSourceAWE/AtmosphericModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/OpenSourceAWE/AtmosphericModels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/OpenSourceAWE/AtmosphericModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/OpenSourceAWE/AtmosphericModels.jl)

## Installation
Install [Julia 1.10](http://www.julialang.org) or later, if you haven't already. You can add AtmosphericModels from  Julia's package manager, by typing 
```julia
using Pkg
pkg"add AtmosphericModels"
``` 
at the Julia prompt.

## Exported types
```julia
AtmosphericModel
@enum ProfileLaw EXP=1 LOG=2 EXPLOG=3 FAST_EXP=4 FAST_LOG=5 FAST_EXPLOG=6
```

## Exported functions
```julia
clear(s::AM)
calc_rho(s::AM, height)
calc_wind_factor(am::AM, height, profile_law::Int64=am.set.profile_law)
```
## Wind profile

<p align="center"><img src="./doc/wind_profile.png" width="500" /></p>

The EXPLOG profile law is the fitted linear combination of the exponential and the log law.

## Usage
```julia
using AtmosphericModels
am = AtmosphericModel()

const profile_law = Int(EXPLOG)
height = 100.0
wf = calc_wind_factor(am, height, profile_law)
```
The result is the factor with which the ground wind speed needs to be multiplied
to get the wind speed at the given height.

## Plot a wind profile
```julia
using AtmosphericModels, Plots
am = AtmosphericModel()

heights = 6:1000
wf = [calc_wind_factor(am, height, Int(EXPLOG)) for height in heights]

plot(heights, wf, legend=false, xlabel="height [m]", ylabel="wind factor")
```

```julia
using AtmosphericModels, ControlPlots
am = AtmosphericModel()
AtmosphericModels.se().alpha = 0.234  # set the exponent of the power law

heights = 6:200
wf = [calc_wind_factor(am, height, Int(EXP)) for height in heights]

plot(heights, wf, xlabel="height [m]", ylabel="wind factor")
```

## Benchmark
```julia
using AtmosphericModels, BenchmarkTools

am = AtmosphericModel()
@benchmark calc_wind_factor(am, height, Int(EXPLOG)) setup=(height=Float64((6.0+rand()*500.0)))
```
|Profile law|time [ns]|
| ---    |:---:|
|EXP     |12   |
|LOG     |16   |
|EXPLOG  |33   |
|FAST_EXP|6.6  |
|FAST_LOG|6.6  |
|FAST_EXPLOG|6.6|

The FAST versions are an approximations with an error of less than $1.5 \cdot 10^{-5}$ and are correct only for the default values of h_ref, z0 and alpha.

## Air density
```julia
using AtmosphericModels, BenchmarkTools
am = AtmosphericModel()
@benchmark calc_rho(am, height) setup=(height=Float64((6.0+rand()*500.0)))
```
This gives 4.85 ns as result. Plot the air density:
```julia
heights = 6:1000
rhos = [calc_rho(am, height) for height in heights]
plot(heights, rhos, legend=false, xlabel="height [m]", ylabel="air density [kg/m³]")
```
<p align="center"><img src="./doc/airdensity.png" width="500" /></p>

## Running the test scripts
First, add TestEnv to your global environment.
```julia
julia
using Pkg
Pkg.add("TestEnv")
exit()
```
Then you can run Julia using this project and run the tests:
```julia
julia --project
using TestEnv
TestEnv.activate()
include("test/bench.jl")
include("calc_approximations.jl")
include("runtests.jl")
```

## Further reading
These models are described in detail in [Dynamic Model of a Pumping Kite Power System](http://arxiv.org/abs/1406.6218).

## Licence
This project is licensed under the MIT License. Please see the below WAIVER in association with the license.

## WAIVER
Technische Universiteit Delft hereby disclaims all copyright interest in the package “AtmosphericModels.jl” (models for airborne wind energy systems) written by the Author(s).

Prof.dr. H.G.C. (Henri) Werij, Dean of Aerospace Engineering

## See also
- [Research Fechner](https://research.tudelft.nl/en/publications/?search=Uwe+Fechner&pageSize=50&ordering=rating&descending=true)
- The application [KiteViewer](https://github.com/ufechner7/KiteViewer)
- the package [KiteUtils](https://github.com/ufechner7/KiteUtils.jl)
- the packages [KiteModels](https://github.com/ufechner7/KiteModels.jl) and [WinchModels](https://github.com/aenarete/WinchModels.jl) and [KitePodModels](https://github.com/aenarete/KitePodModels.jl)
- the packages [KiteControllers](https://github.com/aenarete/KiteControllers.jl) and [KiteViewers](https://github.com/aenarete/KiteViewers.jl)

