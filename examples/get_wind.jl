using AtmosphericModels, KiteUtils

set_data_path("data")
set = load_settings("system.yaml"; relax=true)
am::AtmosphericModel = AtmosphericModel(set)

@info "Ground wind speed: $(am.set.v_wind) m/s"

x, y, z = 20.0, 0.0, 200.0
t = 0.0
vx, vy, vz = get_wind(am, x, y, z, t)
@time get_wind(am, x, y, z, t)
@info "Wind at x=$(x), y=$(y), z=$(z), t=$(t): v_x=$(vx), v_y=$(vy), v_z=$(vz)"
@info "Wind speed: $(sqrt(vx^2 + vy^2 + vz^2)) m/s"

