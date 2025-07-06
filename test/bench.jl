using AtmosphericModels, BenchmarkTools

set_data_path()
set = load_settings("system.yaml"; relax=true)
am::AtmosphericModel = AtmosphericModel(set)

am.set.profile_law=6
@benchmark calc_wind_factor(am, height) setup=(height=Float64((6.0+rand()*500.0)))

# 22ns for 3
#  5ns for 6
