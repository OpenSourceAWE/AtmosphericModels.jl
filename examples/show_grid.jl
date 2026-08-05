using Pkg
if ! ("GLMakie" ∈ keys(Pkg.project().dependencies))
    Pkg.activate("examples")
end
using AtmosphericModels, KiteUtils, GLMakie

set_data_path("data")
set = load_settings("system.yaml"; relax=true)
am::AtmosphericModel = AtmosphericModel(set)

function show_grid(x, y, z)
    """
    x: downwind direction
    z: up
    y: orthogonal to x and z
    """
    fig = Figure(size=(800, 600))
    ax = Axis3(fig[1, 1]; xlabel="X Label", ylabel="Y Label", zlabel="Height [m]",
               title="Show Grid")
    scatter!(ax, vec(x), vec(y), vec(z); markersize=2, color=(:blue, 0.1))
    display(GLMakie.Screen(), fig)
end
ny=50; nx=100; nz=50; z_min=25
am.set.grid= [nx, ny, nz, z_min]
y, x, z = AtmosphericModels.create_grid(am)
show_grid(x, y, z)
