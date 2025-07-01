const GRID_STEP   = 2.0    # Resolution of the grid in x and y direction in meters
const HEIGHT_STEP = 2.0    # Resolution in z direction in meters

function ndgrid(xs, ys, zs)
    X = reshape(xs, :, 1, 1)
    Y = reshape(ys, 1, :, 1)
    Z = reshape(zs, 1, 1, :)
    X = repeat(X, 1, length(ys), length(zs))
    Y = repeat(Y, length(xs), 1, length(zs))
    Z = repeat(Z, length(xs), length(ys), 1)
    return X, Y, Z
end

function create_grid(ny=50, nx=100, nz=50, z_min=25; res=GRID_STEP, height_step=HEIGHT_STEP)
    y_range = range(-ny/2, ny/2, length=Int(ny/res)+1)
    x_range = range(0, nx, length=Int(nx/res)+1)
    z_range = range(z_min, z_min+nz, length=Int(nz/height_step)+1)

    # Create meshgrid (Julia's meshgrid returns in order (x, y, z))
    X, Y, Z = ndgrid(x_range, y_range, z_range)

    return Y, X, Z  # To match the Python (y, x, z) order
end

y, x, z = create_grid(100, 4050, 500, 70)

# Domain dimensions
nx, ny, nz = length(x), length(y), length(z)

# Wave number grid
x_range = range(-nx/2, nx/2-1, length=nx)
y_range = range(-ny/2, ny/2-1, length=ny)
z_range = range(-nz/2, nz/2-1, length=nz)

# m2, m1, m3 = meshgrid(y_range, x_range, z_range)
# println("--> $(size(m1)), $(size(m2)), $(size(m3))")


