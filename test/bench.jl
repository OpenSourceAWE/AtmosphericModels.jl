using AtmosphericModels, BenchmarkTools

set_data_path()
set = load_settings("system.yaml")
am::AtmosphericModel = AtmosphericModel(set)

am.set.profile_law=6
# @benchmark calc_wind_factor(am, height, Val{profile_law}) setup=(height=Float64((6.0+rand()*500.0)))
@benchmark calc_wind_factor(am, height) setup=(height=Float64((6.0+rand()*500.0)))
# @benchmark calc_rho(am, height) setup=(height=Float64((6.0+rand()*500.0)))

