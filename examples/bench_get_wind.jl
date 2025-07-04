using Pkg
if ! ("BenchmarkTools" ∈ keys(Pkg.project().dependencies))
    Pkg.activate("examples")
end
using AtmosphericModels, KiteUtils, BenchmarkTools, Timers

set_data_path("data")
set = load_settings("system.yaml")
am = AtmosphericModel(set)

@info "Ground wind speed: $(am.set.v_wind) m/s"
tic()
wf::WindField = WindField(am, am.set.v_wind)
toc()
x, y, z = 20.0, 25.0, 200.0
t = 0.0
vx, vy, vz = get_wind(wf, am, x, y, z, t)
@btime get_wind(wf, am, x, y, z, t)
@info "Wind speed: $(round(sqrt(vx^2 + vy^2 + vz^2), digits=1)) m/s at $z m height."
# 319.979 ns (11 allocations: 192 bytes) on laptop on battery (without lookup)
# 284.933 ns (19 allocations: 368 bytes) on desktop
# 360.128 ns (19 allocations: 368 bytes) on laptop on battery

