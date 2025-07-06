using AtmosphericModels, KiteUtils

set_data_path("data") 
set = load_settings("system.yaml"; relax=true)
am::AtmosphericModel = AtmosphericModel(set; nowindfield=true)

v_wind_gnd = 5.324

x = range(0, 50, length=25)
y = range(0, 800, length=400)
z = range(0, 200, length=100)

u, v, w = AtmosphericModels.create_windfield(x, y, z; sigma1=1.2)
am = AtmosphericModel(set)
am.set.v_wind = v_wind_gnd
AtmosphericModels.addWindSpeed(am, z, u)