using LinearAlgebra
using FFTW
using KiteUtils
using AtmosphericModels

set_data_path("data")
set = load_settings("system.yaml")
am::AtmosphericModel = AtmosphericModel(set)

const GRID_STEP = 2.0  # Grid resolution in x/y (meters)
const HEIGHT_STEP = 2.0  # Z-resolution (meters)

function create_grid(am, res=GRID_STEP)
    res = Int(res)
    y_range = range(-ny÷2, ny÷2, length=ny÷res + 1)
    x_range = range(0, nx, length=nx÷res + 1)
    z_range = range(z_min, z_min+nz, length=nz÷Int(HEIGHT_STEP) + 1)
    return meshgrid(x_range, y_range, z_range)
end

function meshgrid(x, y, z)
    X = [i for i in x, _ in y, _ in z]
    Y = [j for _ in x, j in y, _ in z]
    Z = [k for _ in x, _ in y, k in z]
    return (X, Y, Z)
end

# Create grid
x, y, z = create_grid(am)

# Get dimensions
nx, ny, nz = size(x, 1), size(y, 2), size(z, 3)

# Wave number discretization
m1_range = range(-nx/2, nx/2 - 1, length=nx)
m2_range = range(-ny/2, ny/2 - 1, length=ny)
m3_range = range(-nz/2, nz/2 - 1, length=nz)

m1, m2, m3 = meshgrid(m1_range, m2_range, m3_range)

# Apply ifftshift with epsilon offset
m1 = ifftshift(m1 .+ 1e-6)
m2 = ifftshift(m2 .+ 1e-6)
m3 = ifftshift(m3 .+ 1e-6)
println("--> $(size(m1)), $(size(m2)), $(size(m3))")