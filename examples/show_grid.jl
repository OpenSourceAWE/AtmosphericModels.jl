using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots, AtmosphericModels

function showGrid(x, y, z)
    """
    x: downwind direction
    z: up
    y: orthogonal to x and z
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(x, y, z; s=2, c="blue", alpha=0.1)
    ax.set_xlabel("X Label")
    ax.set_ylabel("Y Label")
    ax.set_zlabel("Height [m]")
    plt.show()
end

grid = AtmosphericModels.createGrid()
x, y, z = grid
showGrid(x, y, z)
