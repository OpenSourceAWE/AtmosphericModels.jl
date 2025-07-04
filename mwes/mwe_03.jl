using Random
using HypergeometricFunctions:_₂F₁
using LinearAlgebra
using FFTW

const HEIGHT_STEP = 2.0
const GRID_STEP = 2.0

Random.seed!(1234)  # reproducible wind fields

# Hypergeometric function wrapper
# function pfq(z)
#     return real(hyp2f1(1/3, 17/6, 4/3, z))
# end

function pfq(z)
    _₂F₁(1. /3 , 17. /6, 4. /3, z)
end

function ndgrid(xs, ys, zs)
    X = reshape(xs, :, 1, 1)
    Y = reshape(ys, 1, :, 1)
    Z = reshape(zs, 1, 1, :)
    X = repeat(X, 1, length(ys), length(zs))
    Y = repeat(Y, length(xs), 1, length(zs))
    Z = repeat(Z, length(xs), length(ys), 1)
    return X, Y, Z
end

function createWindField(x, y, z; sigma1=nothing, gamma=3.9, ae=0.1, length_scale=33.6)
    # Validate sigma1
    if sigma1 !== nothing
        if !(isa(sigma1, Number) || (isa(sigma1, AbstractVector) && length(sigma1) == 3))
            throw(ArgumentError("The parameter 'sigma1' must be a single value or a 3-component vector."))
        end
    end

    # Check monotonicity
    if x[1] > x[end] || y[1] > y[end] || z[1] > z[end]
        throw(ArgumentError("Values of x, y, and z must be monotonically increasing."))
    end

    # If sigma1 is a scalar, convert to vector for component-wise use
    if sigma1 === nothing
        sigma1_vec = [1.0, 1.0, 1.0]  # default if no sigma1 provided
    elseif isa(sigma1, Number)
        sigma1_vec = [sigma1, sigma1, sigma1]
    else
        sigma1_vec = sigma1
    end

    sigma_iso = 0.55 .* sigma1_vec
    sigma2 = 0.7 .* sigma1_vec
    sigma3 = 0.5 .* sigma1_vec

    nx, ny, nz = size(x)
    # Domain lengths
    Lx = x[end,1,1] - x[1,1,1]
    Ly = y[1,end,1] - y[1,1,1]
    Lz = z[1,1,end] - z[1,1,1]

    # Wave number discretization
    y_range = range(-ny/2, stop=ny/2 - 1, length=ny)
    x_range = range(-nx/2, stop=nx/2 - 1, length=nx)
    z_range = range(-nz/2, stop=nz/2 - 1, length=nz)

    # meshgrid equivalent in Julia: use broadcasting
    m1 = reshape(x_range, nx, 1, 1) .+ 1e-6
    m2 = reshape(y_range, 1, ny, 1) .+ 1e-6
    m3 = reshape(z_range, 1, 1, nz) .+ 1e-6

    # fftshift equivalent: use fftshift from FFTW.jl
    m1 = fftshift(m1, 1)
    m2 = fftshift(m2, 2)
    m3 = fftshift(m3, 3)

    k1 = 2pi .* m1 .* (length_scale / Lx)
    k2 = 2pi .* m2 .* (length_scale / Ly)
    k3 = 2pi .* m3 .* (length_scale / Lz)

    k = sqrt.(k1.^2 .+ k2.^2 .+ k3.^2)

    pfq_term = pfq.(-k.^(-2))
    beta = gamma ./ (k.^(2/3) .* sqrt.(pfq_term))

    k30 = k3 .+ beta .* k1
    k0 = sqrt.(k1.^2 .+ k2.^2 .+ k30.^2)

    E0 = ae * length_scale^(5/3) .* k0.^4 ./ (1 .+ k0.^2).^(17/6)

    # Avoid division by zero by adding a small epsilon
    eps = 1e-14
    denom = k.^2 .* (k1.^2 .+ k2.^2) .+ eps

    C1 = (beta .* k1.^2 .* (k1.^2 .+ k2.^2 .- k3 .* (k3 .+ beta .* k1))) ./ denom
    arctan_arg = (beta .* k1 .* sqrt.(k1.^2 .+ k2.^2))
    arctan_denom = k0.^2 .- (k3 .+ beta .* k1) .* k1 .* beta
    C2 = (k2 .* k0.^2) ./ (k1.^2 .+ k2.^2).^(3/2) .* atan.(arctan_arg, arctan_denom)

    zeta1 = C1 .- k2 ./ k1 .* C2
    zeta2 = C2 .+ k2 ./ k1 .* C1

    B = sigma_iso[1] * sqrt.(2pi^2 * length_scale^3 .* E0 ./ (Lx * Ly * Lz .* k0.^4))

    # Initialize correlation matrix C with dimensions (3,3,nx,ny,nz)
    C = zeros(ComplexF64, 3, 3, nx, ny, nz)

    C[1,1,:,:,:] = B .* k2 .* zeta1
    C[1,2,:,:,:] = B .* (k30 .- k1 .* zeta1)
    C[1,3,:,:,:] = B .* -k2
    C[2,1,:,:,:] = B .* (k2 .* zeta2 .- k30)
    C[2,2,:,:,:] = B .* -k1 .* zeta2
    C[2,3,:,:,:] = B .* k1
    C[3,1,:,:,:] = B .* k2 .* k0.^2 ./ k.^2
    C[3,2,:,:,:] = B .* -k1 .* k0.^2 ./ k.^2
    # C[3,3,:,:,:] remains zero

    # Generate white noise vector n with shape (3, 1, nx, ny, nz)
    n_real = randn(3, 1, nx, ny, nz)
    n_imag = randn(3, 1, nx, ny, nz)
    n = complex.(n_real, n_imag)

    # Compute stochastic field dZ (3, nx, ny, nz)
    dZ = zeros(ComplexF64, 3, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k_ in 1:nz
        # Extract 3x3 matrix C[:,:,i,j,k] and 3x1 vector n[:,1,i,j,k]
        C_mat = reshape(C[:,:,i,j,k_], (3,3))
        n_vec = reshape(n[:,1,i,j,k_], 3)
        dZ[:,i,j,k_] = C_mat * n_vec
    end

    # Inverse FFT and scaling
    # u = nx * ny * nz * real.(ifftn(dZ[1,:,:,:]))
    # v = nx * ny * nz * real.(ifftn(dZ[2,:,:,:]))
    # w = nx * ny * nz * real.(ifftn(dZ[3,:,:,:]))

    u = nx * ny * nz * real.(ifft(dZ[1,:,:,:]))
    v = nx * ny * nz * real.(ifft(dZ[2,:,:,:]))
    w = nx * ny * nz * real.(ifft(dZ[3,:,:,:]))

    if sigma1 !== nothing
        su = std(vec(u))
        sv = std(vec(v))
        sw = std(vec(w))
        u .*= sigma1_vec[1] / su
        v .*= sigma2[2] / sv
        w .*= sigma3[3] / sw
    end

    return u, v, w
end

# function createGrid(ny=50, nx=100, nz=50, z_min=25, res=GRID_STEP)
#     res = Int(res)
#     y_range = range(-ny÷2, stop=ny÷2, length=div(ny,res)+1)
#     x_range = range(0, stop=nx, length=div(nx,res)+1)
#     z_range = range(z_min, stop=z_min+nz, length=div(nz,Int(HEIGHT_STEP))+1)

#     # meshgrid: y, x, z shapes (nx, ny, nz)
#     y = reshape(y_range, 1, length(y_range), 1)
#     x = reshape(x_range, length(x_range), 1, 1)
#     z = reshape(z_range, 1, 1, length(z_range))

#     # Broadcast to get full 3D arrays
#     Y = broadcast(*, ones(length(x_range)), y, ones(length(z_range)))
#     X = broadcast(*, x, ones(length(y_range)), ones(length(z_range)))
#     Z = broadcast(*, ones(length(x_range)), ones(length(y_range)), z)

#     return Y, X, Z
# end
function create_grid(ny=50, nx=100, nz=50, z_min=25; res=GRID_STEP, height_step=HEIGHT_STEP)
    y_range = range(-ny/2, ny/2, length=Int(ny/res)+1)
    x_range = range(0, nx, length=Int(nx/res)+1)
    z_range = range(z_min, z_min+nz, length=Int(nz/height_step)+1)

    # Create meshgrid (Julia's meshgrid returns in order (x, y, z))
    X, Y, Z = ndgrid(x_range, y_range, z_range)

    return Y, X, Z  # To match the Python (y, x, z) order
end

# Example usage
y, x, z = create_grid(10, 20, 10, 5)
u, v, w = createWindField(x, y, z, sigma1=1.0)
nothing
