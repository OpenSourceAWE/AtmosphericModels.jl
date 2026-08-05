using Pkg
if ! ("MakieControlPlots" ∈ keys(Pkg.project().dependencies))
    Pkg.activate("examples")
end

using Statistics, KiteUtils, AtmosphericModels, MakieControlPlots

set_data_path("data")
set = load_settings("system.yaml"; relax=true)
am::AtmosphericModel = AtmosphericModel(set)

function plot_wind_vs_time(am, x=0.0, y=0.0; z=197.3)
    rel_turb = rel_turbo(am)
    @info "Relative turbulence: $rel_turb"
    TIME = Float64[]
    v_wind_x = Float64[]
    v_wind_norm = Float64[]
    for t in range(0.0, stop=600.0, length=600*20)
        push!(TIME, t)
        v_x, v_y, v_z = get_wind(am, x, y, z, t)
        v_wind = sqrt(v_x^2 + v_y^2 + v_z^2)
        push!(v_wind_norm, v_wind)
        push!(v_wind_x, v_x)
        if v_wind < 0.1
            println("Error for x, y, z, t: ", x, y, z, t)
        end
    end
    su = std(v_wind_x)
    mean_val = round(mean(v_wind_x), digits=1)
    turbulence_intensity = round(su / mean_val * 100.0, digits=1)
    println("Mean wind x: $(mean_val) m/s, turbulence intensity: $(turbulence_intensity) %")
    p = plot(TIME, v_wind_x; xlabel="Time [s]", ylabel="Abs. wind speed at $z m height [m/s]",
             fig="Wind speed at I = $(turbulence_intensity) %, z= $(z) m")
    display(p)
end

plot_wind_vs_time(am; z=197.3)
plot_wind_vs_time(am; z=100.0)
nothing

