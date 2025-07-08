using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    Pkg.activate("examples")
end
using AtmosphericModels, KiteUtils, ControlPlots

set_data_path("data")
set = load_settings("system.yaml"; relax=true)
set.v_wind = 3.483*1.085  # wind speed at reference height, average wind speed at Cabauw, NL [m/s]

@info "Ground wind speed: $(set.v_wind) m/s"
am::AtmosphericModel = AtmosphericModel(set)

heights = 10:200
v_w = set.v_wind .* [calc_wind_factor(am, height) for height in heights]

p=plot(heights, v_w,  ylabel="v_wind [m/s]", xlabel="height [m]", fig="Cabauw, NL, EXP law")
display(p)
