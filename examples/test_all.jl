# Expected results for low, medium and high wind speed:
# vw_10m    I 99  I197
# 4.26 m/s  8.5 % 6.3 %
# 6.00 m/s  9.7 % 7.2 %
# 9.20 m/s  9.8 % 7.9 %

using Statistics, KiteUtils, AtmosphericModels, Test

set_data_path("data")
set = load_settings("system.yaml")

I99::Vector{Float64} = [8.5, 9.7, 9.8] # Turbulence intensity at 99 m/s
I197::Vector{Float64} = [6.3, 7.2, 7.9] # Turbulence intensity at 197 m/s

function analyze_windfield(am::AtmosphericModel; z=197.3)
    x = 0.0; y = 0.0
    v_wind_x = Float64[]
    v_wind_norm = Float64[]
    for y in -50:5:50
        for t in range(0.0, stop=600.0, length=601)
            v_x, v_y, v_z = get_wind(am, x, y, z, t)
            v_wind = sqrt(v_x^2 + v_y^2 + v_z^2)
            push!(v_wind_norm, v_wind)
            push!(v_wind_x, v_x)
        end
    end
    su = std(v_wind_norm)
    v_mean = round(mean(v_wind_norm), digits=1)
    turbulence_intensity = round(su / mean(v_wind_norm) * 100.0, digits=2)
    return v_mean, turbulence_intensity 
end

for (i, v_wind_gnd) in pairs(am.set.v_wind_gnds)
    local v_mean, ti, am
    set.v_wind = v_wind_gnd
    am::AtmosphericModel = AtmosphericModel(set)
    v_mean, ti = analyze_windfield(am)
    v_mean_99, ti_99 = analyze_windfield(am; z=99.0)
    @test ti_99 ≈ I99[i]  rtol=0.07
    @test ti    ≈ I197[i] rtol=0.07
    @info "v_mean_99: $v_mean_99 m/s, ti_99: $ti_99 %"
    @info "v_mean_197: $v_mean m/s, ti_197: $ti %"
end