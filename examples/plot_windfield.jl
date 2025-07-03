using Pkg
if ! ("GLMakie" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using GLMakie, AtmosphericModels, KiteUtils

num_points = 100
V_MIN::Float64 = -100.0
V_MAX::Float64 =  100.0
V_WIND::Vector{Float64} = zeros(num_points)

set_data_path("data")
set = load_settings("system.yaml")
am = AtmosphericModel(set)
wf::WindField = WindField(am, am.set.v_wind)
@info "Wind speed at refence height: $(am.set.v_wind) m/s"

function v_wind(pos_x, pos_z, time, num_points)
    pos_y = range(V_MIN, V_MAX, length=num_points)
    for (i, y) in pairs(pos_y)
        vx, vy, vz = get_wind(wf, am, pos_x, y, pos_z, time)
        V_WIND[i] = sqrt(vx^2+vy^2+vz^2)
    end
    return pos_y, V_WIND
end

pos_x = Observable(20.0)   # x-position
pos_z = Observable(200.0) # height
time = Observable(0.0)    # time


function v_wind_y(pos_x, pos_z, time)
    x, y = v_wind(pos_x, pos_z, time, num_points)
    return y
end

y = @lift(v_wind_y($pos_x, $pos_z, $time))
x = range(V_MIN, V_MAX, length=num_points)

# Create the figure and axis
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "pos_y [m]", ylabel = "v_wind [m/s]",
          title = "Wind Field")

# Plot the sine wave
lineplot = lines!(ax, x, y)
ylims!(ax, 0, 20)
xlims!(ax, V_MIN, V_MAX)

# Create sliders for amplitude and frequency
sg = SliderGrid(
    fig[2, 1],
    (label = "pos_x", range = 0:1:500.0, startvalue =  20.0, update_while_dragging=true),
    (label = "pos_z", range = 5:1:500.0, startvalue = 200.0, update_while_dragging=true),
    (label = "time", range = 0.0:1:1000.0, startvalue = 0.0, update_while_dragging=true),
)

# Connect sliders to observables
on(sg.sliders[1].value) do val
    pos_x[] = val
end

on(sg.sliders[2].value) do val
    pos_z[] = val
end

on(sg.sliders[3].value) do val
    time[] = val
end

fig
