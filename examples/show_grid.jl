using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    Pkg.activate("examples")
end
using ControlPlots, AtmosphericModels

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

grid = AtmosphericModels.create_grid()
x, y, z = grid
show_grid(x, y, z)
