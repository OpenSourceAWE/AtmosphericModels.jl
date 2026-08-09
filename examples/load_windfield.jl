using AtmosphericModels, KiteUtils

set_data_path("data") 
set = load_settings("system.yaml"; relax=true)
am::AtmosphericModel = AtmosphericModel(set)

u, v, w, param, v_wind_gnd = AtmosphericModels.load_windfield(am, 5.324)

@info "Windfield dimensions: $(size(u))"
@info "Grid axes: $(AtmosphericModels.grid_axes(am))"
@info "Parameters: $param, generated for $v_wind_gnd m/s"