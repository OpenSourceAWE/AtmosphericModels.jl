using AtmosphericModels, KiteUtils

set_data_path("data") 
set = load_settings("system.yaml"; relax=true)
am::AtmosphericModel = AtmosphericModel(set)

x,y,z,u,v,w,param = AtmosphericModels.load_windfield(am, 5.324)

@info "Windfield dimensions: $(size(x)), $(size(y)), $(size(z))"
@info "Parameters: $param"