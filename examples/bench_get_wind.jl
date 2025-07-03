using Pkg
if ! ("BenchmarkTools" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using AtmosphericModels, KiteUtils, BenchmarkTools

set_data_path("data")
set = load_settings("system.yaml")
am = AtmosphericModel(set)

@info "Ground wind speed: $(am.set.v_wind) m/s"

wf::WindField = WindField(am, am.set.v_wind)
x, y, z = 20.0, 0.0, 200.0
t = 0.0
vx, vy, vz = get_wind(wf, am, x, y, z, t)
@btime get_wind(wf, am, x, y, z, t)
@info "Wind speed: $(round(sqrt(vx^2 + vy^2 + vz^2), digits=1)) m/s at $z m height."

