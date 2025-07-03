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

function turbulence_intensity(vec::Vector)
    TURB=V_WIND_ABS .- mean(V_WIND_ABS)
    v_mean = mean(V_WIND_ABS)
    rms = sqrt(mean(TURB .^ 2))
    I = round(100*rms/v_mean; digits=1)
end
function turbulence_intensity(wf::WindField, height::Number)
    pos_z = height
    pos_x = 0.0
    y = 0.0
    V_ABS = zeros(601)
    for time in 0:600
        vx, vy, vz = get_wind(wf, am, pos_x, y, pos_z, time)
        v_abs = sqrt(vx^2+vy^2+vz^2)
        V_ABS[time+1] = v_abs
    end
    turbulence_intensity(V_ABS)
end
function wind_speed(wf::WindField, height::Number)
    pos_z = height
    pos_x = 0.0
    y = 0.0
    V_ABS = zeros(601)
    for time in 0:600
        vx, vy, vz = get_wind(wf, am, pos_x, y, pos_z, time)
        v_abs = sqrt(vx^2+vy^2+vz^2)
        V_ABS[time+1] = v_abs
    end
    round(mean(V_ABS), digits=1)
end

set_data_path("data")
set = load_settings("system.yaml")
set.v_wind = 10.9
am = AtmosphericModel(set)
wf::WindField = WindField(am, am.set.v_wind)
@info "Wind speed at reference height: $(am.set.v_wind) m/s"


# Create the figure and axis
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "pos_y [m]", ylabel = "v_wind [m/s]",
          title = "Wind Field")

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
    turb=V_WIND_ABS .- mean(V_WIND_ABS)
    v_mean = mean(V_WIND_ABS)
    rms = sqrt(mean(turb .^ 2))
    I = round(100*rms/v_mean; digits=1)
    v_mean= round(v_mean; digits=1)
    ax.title="Wind Field, I = $I %, v_mean = $v_mean m/s"
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
@info "Wind speed at 100m height:           $(wind_speed(wf,100)) m/s"
@info "Wind speed at 200m height:           $(wind_speed(wf,200)) m/s"
@info "Turbulence intensity at 200m height: $(turbulence_intensity(wf, 200)) %"

fig
