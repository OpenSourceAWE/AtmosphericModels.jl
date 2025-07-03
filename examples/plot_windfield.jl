using Pkg
if ! ("GLMakie" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using GLMakie, AtmosphericModels, KiteUtils, Statistics

num_points = 100
V_MIN::Float64 = -100.0
V_MAX::Float64 =  100.0
I::Float64 = 0.0 # turbulence intensity
V_WIND_ABS::Vector{Float64} = zeros(num_points)
V_WIND_X::Vector{Float64} = zeros(num_points)
V_WIND_Y::Vector{Float64} = zeros(num_points)
V_WIND_Z::Vector{Float64} = zeros(num_points)

set_data_path("data")
set = load_settings("system.yaml")
am = AtmosphericModel(set)
wf::WindField = WindField(am, am.set.v_wind)
@info "Wind speed at refence height: $(am.set.v_wind) m/s"

function v_wind(pos_x, pos_z, time, num_points)
    global I
    pos_y = range(V_MIN, V_MAX, length=num_points)
    for (i, y) in pairs(pos_y)
        vx, vy, vz = get_wind(wf, am, pos_x, y, pos_z, time)
        V_WIND_X[i] = vx
        V_WIND_Y[i] = vy
        V_WIND_Z[i] = vz
        V_WIND_ABS[i] = sqrt(vx^2+vy^2+vz^2)
    end
    rms = sqrt(mean(V_WIND_ABS .^ 2))
    I = rms/mean(V_WIND_ABS)
    return V_WIND_ABS
end

pos_x = Observable(20.0)   # x-position
pos_z = Observable(200.0) # height
time = Observable(0.0)    # time


function v_wind_abs(pos_x, pos_z, time)
    v_abs = v_wind(pos_x, pos_z, time, num_points)
    return v_abs
end
function v_wind_x(pos_x, pos_z, time)
    v_wind(pos_x, pos_z, time, num_points)
    return V_WIND_X
end
function v_wind_y(pos_x, pos_z, time)
    v_wind(pos_x, pos_z, time, num_points)
    return V_WIND_Y
end
function v_wind_z(pos_x, pos_z, time)
    v_wind(pos_x, pos_z, time, num_points)
    return V_WIND_Z
end


v_abs = @lift(v_wind_abs($pos_x, $pos_z, $time))
v_x = @lift(v_wind_x($pos_x, $pos_z, $time))
v_y = @lift(v_wind_y($pos_x, $pos_z, $time))
v_z = @lift(v_wind_z($pos_x, $pos_z, $time))
x = range(V_MIN, V_MAX, length=num_points)

# Create the figure and axis
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "pos_y [m]", ylabel = "v_wind [m/s]",
          title = "Wind Field")

# Plot wind speed
line_abs = lines!(ax, x, v_abs)
line_x = lines!(ax, x, v_x)
line_y = lines!(ax, x, v_y)
line_z = lines!(ax, x, v_z)
ylims!(ax, -5, 20)
xlims!(ax, V_MIN, V_MAX)

Legend(fig[1, 2],
    [line_abs, line_x, line_y, line_z],
    ["v_abs", "v_x", "v_y", "v_z"])

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

# rms = sqrt(mean(A .^ 2))

fig
