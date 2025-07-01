using AtmosphericModels

set_data_path("data")
set = load_settings("system.yaml")
am = AtmosphericModel(set)

@info "Ground wind speed: $(am.set.v_wind) m/s"

