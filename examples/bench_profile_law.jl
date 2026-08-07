using Pkg
if ! ("BenchmarkTools" ∈ keys(Pkg.project().dependencies))
    Pkg.activate("examples")
end
using AtmosphericModels, KiteUtils, BenchmarkTools

set_data_path("data")
set = load_settings("system.yaml"; relax=true)
set.heights = [6.0, 20.0, 60.0, 100.0, 200.0]
set.speeds  = [5.324, 6.8, 8.1, 8.9, 9.8]

am::AtmosphericModel = AtmosphericModel(set)
height = 100.0

println("EXP:")
@btime calc_wind_factor(am, height, Int(EXP))
println("LOG:")
@btime calc_wind_factor(am, height, Int(LOG))
println("EXPLOG:")
@btime calc_wind_factor(am, height, Int(EXPLOG))
println("CUSTOM_EXP:")
@btime calc_wind_factor(am, height, Int(CUSTOM_EXP))
println("CUSTOM_LOG:")
@btime calc_wind_factor(am, height, Int(CUSTOM_LOG))

# CUSTOM_JET fits 5 coefficients, so it needs jet-shaped data to converge realistically.
c_bg, a_bg, U_J, z_c, sigma = 3.0, 0.15, 4.0, 200.0, 40.0
am.set.heights = [10.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 600.0]
am.set.speeds  = [c_bg * z^a_bg + U_J * exp(-(z - z_c)^2 / (2 * sigma^2)) for z in am.set.heights]

# am caches the last fit in am.jet_cache, reused while am.set.heights/speeds are unchanged.
println("CUSTOM_JET (cache miss, fresh am per call):")
@btime calc_wind_factor(am2, height, Int(CUSTOM_JET)) setup=(am2 = AtmosphericModel(am.set; nowindfield=true)) evals=1
println("CUSTOM_JET (cached):")
@btime calc_wind_factor(am, height, Int(CUSTOM_JET))
# on desktop:
# EXP:                             28.362 ns (1 allocation: 16 bytes)
# LOG:                             24.769 ns (1 allocation: 16 bytes)
# EXPLOG:                          34.099 ns (1 allocation: 16 bytes)
# CUSTOM_EXP:                     126.950 ns (11 allocations: 432 bytes)
# CUSTOM_LOG:                      80.052 ns (9 allocations: 336 bytes)
# CUSTOM_JET (cache miss):        9.500 μs (287 allocations: 24.73 KiB)
# CUSTOM_JET (cached):            34.693 ns (0 allocations: 0 bytes)
