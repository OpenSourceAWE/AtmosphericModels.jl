# AtmosphericModels

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://OpenSourceAWE.github.io/AtmosphericModels.jl/dev)
[![Build Status](https://github.com/OpenSourceAWE/AtmosphericModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/OpenSourceAWE/AtmosphericModels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/OpenSourceAWE/AtmosphericModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/OpenSourceAWE/AtmosphericModels.jl)

## Installation
Install [Julia 1.10](https://ufechner7.github.io/2024/08/09/installing-julia-with-juliaup.html) or later, if you haven't already. You can add AtmosphericModels from  Julia's package manager, by typing 

First, create a new Julia project:
```bash
mkdir test
cd test
julia --project=.
```

You can now add AtmosphericModels from  Julia's package manager, by typing 
```julia
using Pkg
pkg"add AtmosphericModels"
``` 
at the Julia prompt.

## Running the tests
Launch Julia using this project and run the tests:
```julia
julia --project
julia> using Pkg
julia> Pkg.test("AtmosphericModels")
```

## Running the examples
If you check out the project using git, you can more easily run the examples:
```bash
git clone https://github.com/OpenSourceAWE/AtmosphericModels.jl
cd AtmosphericModels.jl
```
Launch Julia using this project with `julia --project` and run the example menu:
```julia
include("examples/menu.jl")
```
The first time will take some time, because the graphic libraries will get installed, the second time it is fast.

## Usage
### Calculate the height dependant wind speed
Make sure that the folder `data` exist and contains the files `system_nearshore.yaml` and `settings_nearshore.yaml`.
These configuration files contain the wind profile parameters, fitted to the near shore location Maasvlakte, NL
on a specific day.

```julia
using AtmosphericModels, KiteUtils
set_data_path("data")
set = load_settings("system.yaml"; relax=true)
am = AtmosphericModel(set)

height = 100.0
wf = calc_wind_factor(am, height)
```
The result is the factor with which the ground wind speed needs to be multiplied
to get the wind speed at the given height.

## Using the turbulent wind field
You can get a wind vector as function of x,y,z and time using the following code:
```julia
using AtmosphericModels, KiteUtils

set_data_path("data")
set = load_settings("system.yaml"; relax=true)
am::AtmosphericModel = AtmosphericModel(set)

@info "Ground wind speed: $(am.set.v_wind) m/s"

wf::WindField = WindField(am, am.set.v_wind)
x, y, z = 20.0, 0.0, 200.0
t = 0.0
vx, vy, vz = get_wind(wf, am, x, y, z, t)
@time get_wind(am, x, y, z, t)
@info "Wind at x=$(x), y=$(y), z=$(z), t=$(t): v_x=$(vx), v_y=$(vy), v_z=$(vz)"
@info "Wind speed: $(sqrt(vx^2 + vy^2 + vz^2)) m/s"
```
It is suggested to check out the code using git before executing this example,
because it requires that a data directory with the correct files `system.yaml`
and `settings.yaml` exists. See below how to do that.

## Plot a wind profile
```julia
using AtmosphericModels, KiteUtils, ControlPlots
am = AtmosphericModel(se())

heights = 6:1000
wf = [calc_wind_factor(am, height, Int(EXPLOG)) for height in heights]

plot(heights, wf, xlabel="height [m]", ylabel="wind factor", fig="Nearshore")
```
![Wind profile nearshore](docs/src/nearshore.png)
```julia
using AtmosphericModels, ControlPlots, KiteUtils
am = AtmosphericModel(se())
AtmosphericModels.se().alpha = 0.234  # set the exponent of the power law

heights = 6:200
wf = [calc_wind_factor(am, height, Int(EXP)) for height in heights]

plot(heights, wf, xlabel="height [m]", ylabel="wind factor", fig="Onshore")
```

## Air density
```julia
using AtmosphericModels, BenchmarkTools, KiteUtils
am = AtmosphericModel(se())
@benchmark calc_rho(am, height) setup=(height=Float64((6.0+rand()*500.0)))
```
This gives 4.85 ns as result. Plot the air density:
```julia
heights = 6:1000
rhos = [calc_rho(am, height) for height in heights]
plot(heights, rhos, legend=false, xlabel="height [m]", ylabel="air density [kg/m³]")
```
<p align="center"><img src="./docs/src/airdensity.png" width="500" /></p>

## Further reading
Please read the [documentation](https://OpenSourceAWE.github.io/AtmosphericModels.jl/dev). At the end of the documentation ([here](https://opensourceawe.github.io/AtmosphericModels.jl/dev/wind_field/#References)) you find references to the scientific literature.

## License
This project is licensed under the MIT License. Please see the below `Copyright notice` in association with the licenses that can be found in the file [LICENSE](LICENSE) in this folder.

## Copyright notice
Technische Universiteit Delft hereby disclaims all copyright interest in the package “AtmosphericModels.jl” (models for airborne wind energy systems) written by the Author(s).

Prof.dr. H.G.C. (Henri) Werij, Dean of Aerospace Engineering, Technische Universiteit Delft.

See the copyright notices in the source files, and the list of authors in [AUTHORS.md](AUTHORS.md).

## See also
- [Research Fechner](https://research.tudelft.nl/en/publications/?search=Uwe+Fechner&pageSize=50&ordering=rating&descending=true)
- The application [KiteViewer](https://github.com/ufechner7/KiteViewer)
- the package [KiteUtils](https://github.com/ufechner7/KiteUtils.jl)
- the packages [KiteModels](https://github.com/ufechner7/KiteModels.jl) and [WinchModels](https://github.com/aenarete/WinchModels.jl) and [KitePodModels](https://github.com/aenarete/KitePodModels.jl)
- the packages [KiteControllers](https://github.com/aenarete/KiteControllers.jl) and [KiteViewers](https://github.com/aenarete/KiteViewers.jl)

