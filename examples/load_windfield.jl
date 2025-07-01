using AtmosphericModels, KiteUtils

set_data_path("data") 
set = load_settings("system.yaml")
am = AtmosphericModel(set)

x,y,z,u,v,w,param = AtmosphericModels.load_windfield(am, 5.324)

println("Windfield dimensions: $(size(x)), $(size(y)), $(size(z))")