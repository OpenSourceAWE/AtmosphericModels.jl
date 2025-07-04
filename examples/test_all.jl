# test the wind field for all three test cases

using ControlPlots, Statistics, KiteUtils, AtmosphericModels

set_data_path("data")
set = load_settings("system.yaml")
am = AtmosphericModel(set)

function analyze_windfield(wf::WindField, am; z=197.3)
    x = 0.0; y = 0.0
    v_wind_x = Float64[]
    v_wind_norm = Float64[]
    for t in range(0.0, stop=600.0, length=600*20)
        v_x, v_y, v_z = get_wind(wf, am, x, y, z, t)
        v_wind = sqrt(v_x^2 + v_y^2 + v_z^2)
        push!(v_wind_norm, v_wind)
        push!(v_wind_x, v_x)
        if v_wind < 0.1
            println("Error for x, y, z, t: ", x, y, z, t)
        end
    end
    su = std(v_wind_x)
    v_mean = round(mean(v_wind_x), digits=1)
    turbulence_intensity = round(su / mean(v_wind_x) * 100.0, digits=1)
    return v_mean, turbulence_intensity 
end

for v_wind_gnd in am.set.v_wind_gnds
    local v_mean, ti, wf
    am.set.v_wind = v_wind_gnd
    wf = WindField(am, am.set.v_wind)
    v_mean, ti = analyze_windfield(wf, am)
    v_mean_100, ti_100 = analyze_windfield(wf, am; z=100.0)
    @info "v_mean_100: $v_mean_100 m/s, ti_100: $ti %"
    @info "v_mean_197: $v_mean m/s, ti_197: $ti %"
end