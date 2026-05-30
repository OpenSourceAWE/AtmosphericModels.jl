using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    Pkg.activate("examples")
end
using AtmosphericModels, KiteUtils, ControlPlots

set_data_path("data")
set = load_settings("system.yaml"; relax=true)
set.profile_law = 0

@info "Ground wind speed: $(set.v_wind) m/s"
am::AtmosphericModel = AtmosphericModel(set=set)

heights = 6:1000
v_w = set.v_wind .* [calc_wind_factor(am, height) for height in heights]

p=ControlPlots.plot(v_w, heights, xlabel="v_wind [m/s]", ylabel="height [m]", fig="No wind shear")
display(p)
