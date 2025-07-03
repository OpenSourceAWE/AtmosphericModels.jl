using Pkg
if ! ("GLMakie" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using GLMakie

function figure_eight_path(A, B, C, D, x0, y0, theta, num_points)
    t = range(0, 2π, length=num_points)
    x = A * sin.(t)
    y = B * sin.(t) .* cos.(t) .+ C .* cos.(t) .+ D .* cos.(2t)
    y = y / (maximum(y)/(B/2))
    # Apply rotation
    x_rot = x .* cos(theta) .- y .* sin(theta) 
    y_rot = x .* sin(theta) .+ y .* cos(theta)
    # Apply translation
    x_final = x_rot .+ x0
    y_final = y_rot .+ y0
    return x_final, y_final
end

A = Observable(20.0)      # width of the figure-eight
B = Observable(10.0)      # height of the figure-eight
C = Observable(0.0)       # size of right part of the figure-eight
D = Observable(0.0)       # asymmetry factor
x0 = Observable(0.0)      # center x-coordinate
y0 = Observable(25.0)     # center y-coordinate
theta = Observable(0.0)     # rotation angle in degrees
num_points = 200

function figure_eight_y(A, B, C, D, x0, y0, theta)
    x, y = figure_eight_path(A, B, C, D, x0, y0, deg2rad(theta), num_points)
    return y
end

function figure_eight_x(A, B, C, D, x0, y0, theta)
    x, y = figure_eight_path(A, B, C, D, x0, y0, deg2rad(theta), num_points)
    return x
end

y = @lift(figure_eight_y($A, $B, $C, $D, $x0, $y0, $theta))
x = @lift(figure_eight_x($A, $B, $C, $D, $x0, $y0, $theta))

# Create the figure and axis
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "azimuth [°]", ylabel = "elevation [°]",
          title = "Figure of Eight Path")

# Plot the sine wave
lineplot = lines!(ax, x, y)
ylims!(ax, 10, 40)
xlims!(ax, -30, 30)

# Create sliders for amplitude and frequency
sg = SliderGrid(
    fig[2, 1],
    (label = "A (width)", range = 10:0.01:30.0, startvalue = 20.0),
    (label = "B (height)", range = 5:0.01:20.0, startvalue = 10.0),
    (label = "C (right_size)", range = -2.0:0.01:2.0, startvalue = 0.0),
    (label = "D (asymmetry)", range = -3:0.01:3.0, startvalue = 0.0),
    (label = "x0", range = -10:0.01:10.0, startvalue = 0.0),
    (label = "y0", range = 20:0.01:30.0, startvalue = 25.0),
    (label = "theta [°]", range = -90:0.01:90, startvalue = 0.0)
)

# Connect sliders to observables
on(sg.sliders[1].value) do val
    A[] = val
end

on(sg.sliders[2].value) do val
    B[] = val
end

on(sg.sliders[3].value) do val
    C[] = val
end

on(sg.sliders[4].value) do val
    D[] = val
end

on(sg.sliders[5].value) do val
    x0[] = val
end

on(sg.sliders[6].value) do val
    y0[] = val
end

on(sg.sliders[7].value) do val
    theta[] = val
end



fig
