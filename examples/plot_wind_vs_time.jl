using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    Pkg.activate("examples")
end

using ControlPlots, Statistics, KiteUtils, AtmosphericModels

const REL_TURB = [0.342, 0.465, 0.583]
const TEST = 2  # Index for the relative turbulence to use

set_data_path("data")
set = load_settings("system.yaml")
am = AtmosphericModel(set)

@info "Ground wind speed: $(am.set.v_wind) m/s"

wf::WindField = WindField(am, am.set.v_wind)
x, y, z = 20.0, 0.0, 200.0
t = 0.0

function plot_wind_vs_time(wf::WindField, am, x=0.0, y=0.0, z=197.3, rel_turb=REL_TURB[TEST])
    println("Relative turbulence: ", rel_turb)
    fig = plt.figure()
    TIME = Float64[]
    v_wind_x = Float64[]
    v_wind_norm = Float64[]
    for t in range(0.0, stop=600.0, length=600*20)
        push!(TIME, t)
        v_x, v_y, v_z = get_wind(wf, am, x, y, z, t, rel_turb=rel_turb)
        v_wind = sqrt(v_x^2 + v_y^2 + v_z^2)
        push!(v_wind_norm, v_wind)
        push!(v_wind_x, v_x)
        if v_wind < 0.1
            println("Error for x, y, z, t: ", x, y, z, t)
        end
    end
    su = std(v_wind_x)
    mean_val = mean(v_wind_x)
    println("Mean wind x, standard deviation, turbulence intensity [%]: ", mean_val, ", ", su, ", ", su/mean_val * 100.0)
    plt.plot(TIME, v_wind_x, label = "Abs. wind speed at 197.3 m [m/s]", color="black")
    plt.xlabel("Time [s]")
    plt.ylabel("Abs. wind speed at 197.3 m height [m/s]")
    plt.legend(loc="upper right")
end

plot_wind_vs_time(wf, am)
nothing

