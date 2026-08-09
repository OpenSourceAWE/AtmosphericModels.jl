using Pkg
if ! ("BenchmarkTools" ∈ keys(Pkg.project().dependencies))
    Pkg.activate("examples")
end
using AtmosphericModels, KiteUtils, BenchmarkTools, Timers

set_data_path("data")
set = load_settings("system.yaml"; relax=true)

@info "Ground wind speed: $(set.v_wind) m/s"
tic()
am::AtmosphericModel = AtmosphericModel(set)
toc()
x, y, z = 20.0, 25.0, 200.0
t = 0.0
vx, vy, vz = get_wind(am, x, y, z, t)
@btime get_wind(am, x, y, z, t)
@info "Wind speed: $(round(sqrt(vx^2 + vy^2 + vz^2), digits=1)) m/s at $z m height."
# 319.979 ns (11 allocations: 192 bytes) on laptop on battery (without lookup)
# 284.933 ns (19 allocations: 368 bytes) on desktop
# 360.128 ns (19 allocations: 368 bytes) on laptop on battery

# The vector method computes the position-independent part (rel_turbo, the wind direction,
# the settings lookups) once per call instead of once per position.
positions = [KiteUtils.SVec3(x + 3i, y - 2i, z + i) for i in 1:100]
res = Vector{KiteUtils.SVec3}(undef, length(positions))
@assert get_wind(am, positions, t) == [KiteUtils.SVec3(get_wind(am, p[1], p[2], p[3], t)) for p in positions]
println("scalar get_wind in a loop, 100 positions:")
@btime [KiteUtils.SVec3(get_wind($am, p[1], p[2], p[3], $t)) for p in $positions]
println("get_wind(am, positions, t):")
@btime get_wind($am, $positions, $t)
println("get_wind!(res, am, positions, t):")
@btime get_wind!($res, $am, $positions, $t)
#   7.214 µs (203 allocations: 10.26 KiB) loop      on laptop on battery
#   3.073 µs (5 allocations: 2.52 KiB)    vector    on laptop on battery
#   3.023 µs (2 allocations: 80 bytes)    in-place  on laptop on battery

