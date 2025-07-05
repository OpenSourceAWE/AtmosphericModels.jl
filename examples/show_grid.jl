using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    Pkg.activate("examples")
end
using ControlPlots, AtmosphericModels, KiteUtils

set_data_path("data")
set = load_settings("system.yaml")
am::AtmosphericModel = AtmosphericModel(set)

function show_grid(x, y, z)
    """
    x: downwind direction
    z: up
    y: orthogonal to x and z
    """
    fig = plt.figure("Show Grid")
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(x, y, z; s=2, c="blue", alpha=0.1)
    ax.set_xlabel("X Label")
    ax.set_ylabel("Y Label")
    ax.set_zlabel("Height [m]")
    plt.show()
end

y, x, z = AtmosphericModels.create_grid(am)
show_grid(x, y, z)
