using AtmosphericModels, KiteUtils

set_data_path("data") 
set = load_settings("system.yaml")
am = AtmosphericModel(set)

v_wind_gnd = 5.324

am = AtmosphericModel(set)
am.set.v_wind = v_wind_gnd

new_windfield(am, v_wind_gnd)
nothing
