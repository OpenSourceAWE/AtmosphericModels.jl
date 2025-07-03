using Pkg
if ! ("GLMakie" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using GLMakie

function v_wind(pos_x, pos_z, time, num_points)
    x = range(-200.0,200.0, length=num_points)
    y = 10*sin.(0.1x)
    return x, y
end

pos_x = Observable(20.0)   # x-position
pos_z = Observable(200.0) # height
time = Observable(0.0)    # time
num_points = 200

function v_wind_y(pos_x, pos_z, time)
    x, y = v_wind(pos_x, pos_z, time, num_points)
    return y
end

function v_wind_x(pos_x, pos_z, time)
    x, y = v_wind(pos_x, pos_z, time, num_points)
    return x
end

y = @lift(v_wind_y($pos_x, $pos_z, $time))
x = @lift(v_wind_x($pos_x, $pos_z, $time))

# Create the figure and axis
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "pos_y [m]", ylabel = "v_wind [m/s]",
          title = "Wind Field")

# Plot the sine wave
lineplot = lines!(ax, x, y)
ylims!(ax, 0, 15)
xlims!(ax, -200, 200)

# Create sliders for amplitude and frequency
sg = SliderGrid(
    fig[2, 1],
    (label = "pos_x", range = 0:0.01:500.0, startvalue =  20.0),
    (label = "pos_z", range = 5:0.01:500.0, startvalue = 200.0),
    (label = "time", range = 0.0:0.01:1000.0, startvalue = 0.0),
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
