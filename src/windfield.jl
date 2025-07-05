"""
This module provides methods to create a turbulent wind field and to read 
the actual wind velocity vector as function of the 3d position
vector and of the time.
It additionally provides functions for plotting the rel. turbulence as function
of the height.
This code is based on the Matlab module Mann.m written by René Bos.
The code is based on the following papers:
 - Mann, Jakob (1994). The Spatial Structure of Neutral Atmospheric
   Surface-Layer Turbulence. Journal of Fluid Mechanics 273,
   pp. 141-168.
 - Mann, Jakob. (1998). Wind Field Simulation. Probabilistic
   Engineering Mechanics 13(4), pp. 269-282.
"""

# TODO: the following values are hardcoded, but should be set in the settings
const V_WIND_GND  = 8.0    # Default value, change as needed
const GRID_STEP   = 2.0    # Resolution of the grid in x and y direction in meters
const HEIGHT_STEP = 2.0    # Resolution in z direction in meters
const SRL = StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}

Base.@kwdef struct WindField
    x_max::Float64 = NaN
    x_min::Float64 = NaN
    y_max::Float64 = NaN
    y_min::Float64 = NaN
    z_max::Float64 = NaN
    z_min::Float64 = NaN
    last_speed::Float64 = 0.0
    valid::Bool = false
    x::Union{SRL, Array{Float64, 3}}
    y::Union{SRL, Array{Float64, 3}}
    z::Union{SRL, Array{Float64, 3}}
    u::Array{Float64, 3}
    v::Array{Float64, 3}
    w::Array{Float64, 3}
    param::Vector{Float64} = [0, 0] # [alpha, v_wind_gnd]
end
function Base.getproperty(wf::WindField, sym::Symbol)
    if sym == :y_range
        getproperty(wf, :y_max) - getproperty(wf, :y_min) 
    elseif sym == :x_range
        getproperty(wf, :x_max) - getproperty(wf, :x_min)
    else
        getfield(wf, sym)
    end
end
function WindField(am, speed; prn=true)
    try
        last_speed = 0.0
        prn && @info "Loading wind field... $speed m/s"
        x, y, z, u, v, w, param = load_windfield(am, speed)
        valid = true
        x_max = maximum(x)
        x_min = minimum(x)
        y_max = maximum(y)
        y_min = minimum(y)
        z_max = maximum(z)
        z_min = minimum(z)
        return WindField(x_max, x_min, y_max, y_min, z_max, z_min, last_speed, valid, x, y, z, u, v, w, param)
    catch
        @error "Error reading wind field!"
        return nothing
    end
end

function pfq(z)
    _₂F₁(1. /3 , 17. /6, 4. /3, z)
end

function calc_sigma1(am, v_wind_gnd)
    v_height = v_wind_gnd * calc_wind_factor(am, am.set.avg_height, Val{Int(EXP)}) 
    am.set.i_ref * (0.75 * v_height + 5.6)
end

function rel_turbo(am::AtmosphericModel, v_wind = am.set.v_wind)
    # Find the closest relative turbulence value for a given ground wind speed
    min_dist, idx = findmin(abs.(am.set.v_wind_gnds .- v_wind))
    return am.set.rel_turbs[idx]
end

"""
Find 2^n that is equal to or greater than i.
"""
function nextpow2(i)
    n = 1
    while n < i
        n *= 2
    end
    n
end

function calc_full_name(v_wind_gnd; basename="windfield_4050_500", rel_sigma=1.0)
    path = get_data_path() * "/"
    name = basename * "_" * @sprintf("%.1f", rel_sigma)
    name *= "_" * @sprintf("%.1f", v_wind_gnd)
    return path * name
end

function save(am, x, y, z, u, v, w, param; basename="windfield_4050_500", v_wind_gnd=V_WIND_GND)
    fullname = calc_full_name(v_wind_gnd; basename, rel_sigma=am.set.use_turbulence)
    @info "Saving wind field to: $fullname.npz"
    # Save as compressed .npz
    NPZ.npzwrite(fullname * ".npz", Dict(
        "x" => x,
        "y" => y,
        "z" => z,
        "u" => u,
        "v" => v,
        "w" => w,
        "param" => param
    ))
end

function load(; basename="windfield_4050_500", v_wind_gnd=8.0)
    fullname = calc_full_name(v_wind_gnd, basename=basename)
    npzfile = NPZ.npzread(fullname * ".npz")
    return npzfile["x"], npzfile["y"], npzfile["z"], npzfile["u"], npzfile["v"], npzfile["w"], npzfile["param"]
end

function load_windfield(am::AtmosphericModel, speed)
    # Find the index of the closest wind speed
    idx = findmin(abs.(am.set.v_wind_gnds .- speed))[2]
    return load(v_wind_gnd = am.set.v_wind_gnds[idx])
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

"""
    create_grid(ny=50, nx=100, nz=50, z_min=25; res=GRID_STEP, height_step=HEIGHT_STEP)

Creates a 3D grid for the wind field model.

# Arguments
- `ny=50`:    Number of grid points in the y-direction.
- `nx=100`:   Number of grid points in the x-direction.
- `nz=50`:    Number of grid points in the z-direction (vertical).
- `z_min=25`: Minimum height (starting z value) of the grid.

# Keyword Arguments
- `res=GRID_STEP`:           Horizontal grid resolution (distance between points in x and y)   [m]
- `height_step=HEIGHT_STEP`: Vertical grid resolution (distance between points in z) [m]

# Returns
A data structure representing the generated 3D grid.
"""
function create_grid(ny=50, nx=100, nz=50, z_min=25; res=GRID_STEP, height_step=HEIGHT_STEP)
    y_range = range(-ny/2, ny/2, length=Int(ny/res)+1)
    x_range = range(0, nx, length=Int(nx/res)+1)
    z_range = range(z_min, z_min+nz, length=Int(nz/height_step)+1)

    # Create meshgrid (Julia's meshgrid returns in order (x, y, z))
    X, Y, Z = ndgrid(x_range, y_range, z_range)

    return Y, X, Z  # To match the Python (y, x, z) order
end

function meshgrid(x, y, z)
    X = [i for i in x, _ in y, _ in z]
    Y = [j for _ in x, j in y, _ in z]
    Z = [k for _ in x, _ in y, k in z]
    return (X, Y, Z)
end

function create_windfield(x, y, z; sigma1=nothing, gamma=3.9, ae=0.1, length_scale=33.6)
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

"""
    create_windfield(x::AbstractVector, y::AbstractVector, z::AbstractVector;
                        sigma1::Union{Nothing, Real, AbstractVector}=nothing,
                        gamma::Real=3.9, ae::Real=0.1, length_scale::Real=33.6)

## Parameters:
- sigma1: Target value(s) for the turbulence intensity. This can either be a single value, in
          which case the statistics are corrected based on the longitudinal component only, or a vector where
          each u,v,w-component can be defined separately.
- gamma:  wind shear (zero for isotropic turbulence, 3.9 for IEC wind profile)
- ae:     Coefficient of the inertial cascade, ae = a*e^(2/3),
          where a = 1.7 is the three-dimensional Kolmogorov
          constant and e = the mean TKE dissipation rate [m/s].

Performance: Python:  17 seconds for a field of 50x800x200 m with 2 m resolution
             Matlab:  17 seconds
"""
function create_windfield_(x::AbstractArray, y::AbstractArray, z::AbstractArray;
                        sigma1::Union{Nothing, Real, AbstractVector}=nothing,
                        gamma::Real=3.9, ae::Real=0.1, length_scale::Real=33.6)
    # Validate inputs
    if sigma1 !== nothing
        if !(sigma1 isa Real) && length(sigma1) ≠ 3
            throw(ArgumentError("sigma1 must be scalar or 3-element vector"))
        end
    end
       
    # Standard deviations
    sigma_iso = 0.55 * sigma1
    sigma2 = 0.7 * sigma1
    sigma3 = 0.5 * sigma1
    
    # Domain dimensions
    nx, ny, nz = size(x, 1), size(y, 2), size(z, 3)
    Lx, Ly, Lz = x[end] - x[1], y[end] - y[1], z[end] - z[1]

    # Wave number discretization
    m1_range = range(-nx/2, nx/2 - 1, length=nx)
    m2_range = range(-ny/2, ny/2 - 1, length=ny)
    m3_range = range(-nz/2, nz/2 - 1, length=nz)
    m1, m2, m3 = meshgrid(m1_range, m2_range, m3_range)
    
    # Apply ifftshift with epsilon offset
    m1 = ifftshift(m1 .+ 1e-6)
    m2 = ifftshift(m2 .+ 1e-6)
    m3 = ifftshift(m3 .+ 1e-6)

     # Wave number vectors
    k1 = 2π * m1 * (length_scale / Lx)
    k2 = 2π * m2 * (length_scale / Ly)
    k3 = 2π * m3 * (length_scale / Lz)

    # Wave number magnitude
    k = sqrt.(k1.^2 .+ k2.^2 .+ k3.^2)
            
    # Non-dimensional distortion time
    pfq_term = pfq.(-k.^-2)  # Assumes pfq is defined elsewhere
    β = @. gamma / (k^(2/3) * sqrt(pfq_term))
    
    # Initial wave vectors
    k30 = @. k3 + β * k1
    k0 = @. sqrt(k1^2 + k2^2 + k30^2)
    
    # Energy spectrum
    E0 = @. ae * length_scale^(5/3) * k0^4 / (1 + k0^2)^(17/6)
    
    # Correlation matrix components
    C1 = @. (β * k1^2 * (k1^2 + k2^2 - k3 * (k3 + β * k1))) / (k^2 * (k1^2 + k2^2))
    C2 = @. (k2 * k0^2) / (k1^2 + k2^2)^(3/2) * atan((β * k1 * sqrt(k1^2 + k2^2)), (k0^2 - (k3 + β * k1) * k1 * β))
    
    ζ1 = @. C1 - (k2 / k1) * C2
    ζ2 = @. C2 + (k2 / k1) * C1
    
    # Amplitude factor
    B = @. sigma_iso * sqrt(2π^2 * length_scale^3 * E0 / (Lx * Ly * Lz * k0^4))
    
    # Correlation tensor
    C = Array{ComplexF64}(undef, 3, 3, nx, ny, nz)
    C[1,1,:,:,:] = @. B * k2 * ζ1
    C[1,2,:,:,:] = @. B * (k30 - k1 * ζ1)
    C[1,3,:,:,:] = @. -B * k2
    C[2,1,:,:,:] = @. B * (k2 * ζ2 - k30)
    C[2,2,:,:,:] = @. -B * k1 * ζ2
    C[2,3,:,:,:] = @. B * k1
    C[3,1,:,:,:] = @. B * k2 * k0^2 / k^2
    C[3,2,:,:,:] = @. -B * k1 * k0^2 / k^2
    
    # White noise field
    n_real = randn(ComplexF64, 3, 1, nx, ny, nz)
    n_imag = randn(ComplexF64, 3, 1, nx, ny, nz)
    n = n_real + im * n_imag
    
    # Stochastic field (vectorized)
    dZ = similar(n, (3, nx, ny, nz))
    @inbounds for i in 1:nx, j in 1:ny, k in 1:nz
        C_slice = @view C[:, :, i, j, k]
        n_slice = @view n[:, :, i, j, k]
        dZ[:, i, j, k] = C_slice * n_slice
    end
    
    # Inverse FFT and scaling
    u = real(ifft(dZ[1,:,:,:])) * nx * ny * nz
    v = real(ifft(dZ[2,:,:,:])) * nx * ny * nz
    w = real(ifft(dZ[3,:,:,:])) * nx * ny * nz
    
    # Normalize variances if requested
    if sigma1 !== nothing
        u .*= sigma1 / std(u)
        v .*= sigma2 / std(v)
        w .*= sigma3 / std(w)
    end
    
    return u, v, w
end

function addWindSpeed(am::AtmosphericModel, z, u)
    """
    Modify the velocity component u such that the average wind speed, calculated according
    to the given wind profile, is added.
    """
    for i in axes(z, 3)
        height = z[1, 1, i]
        v_wind = am.set.v_wind * calc_wind_factor(am, height, am.set.profile_law)
        u[:, :, i] .+= v_wind
    end
end

function load(wf::WindField, speed)
    global ALPHA, V_WIND_GND
    if speed == wf.last_speed
        return
    end
    println("Loading wind field ... $speed m/s")
    wf.last_speed = speed
    wf.x, wf.y, wf.z, wf.u, wf.v, wf.w, wf.param = load_windfield(speed)
    wf.valid = true
    ALPHA = wf.param[0]
    V_WIND_GND = wf.param[1]
    #  self.u_pre, self.v_pre, self.w_pre = [ndimage.spline_filter(item, order=3) for item \
    #                                                                             in [self.u, self.v, self.w]]
    nothing
end

function get_wind(wf::WindField, am::AtmosphericModel, x, y, z, t; interpolate=false)
    """ 
    Return the wind vector for a given position and time. Linear interpolation in x, y and z.
    """
    @assert z >= 5.0 "Height must be at least 5 m"
    if z < 10.0
        z = 10.0
    end
    @assert t >= 0.0 "Time must be non-negative"
    rel_turb = rel_turbo(am)  
    
    # duplicate the wind field in x and y direction
    while x < wf.x_min
        x += wf.x_range
    end
    while y > wf.y_max
        y -= wf.y_range
    end
    while y < wf.y_min
        y += wf.y_range
    end
    
    y1 = ((y + wf.y_range / 2) / GRID_STEP)
    v_wind_height = am.set.v_wind * calc_wind_factor(am, z, am.set.profile_law)
    # v_wind_height = 11.0
    
    x1 = (x + t * v_wind_height) / GRID_STEP
    while x1 > size(wf.u, 1) - 1
        x1 -= size(wf.u, 1) - 1
    end
    x1 = Int(round(x1))+1
    y1 = Int(round(y1))+1 
    
    z1 = z / HEIGHT_STEP
    if z1 > size(wf.u, 3) - 1
        z1 = size(wf.u, 3) - 1
    elseif z1 < 0
        z1 = 0
    end
    z1 = Int(round(z1))+1
    
    if interpolate
        # x_wind = ndimage.map_coordinates(wf.u, [[x1], [y1], [z1]], order=3, prefilter=false)
        # y_wind = ndimage.map_coordinates(wf.v, [[x1], [y1], [z1]], order=3, prefilter=false)
        # z_wind = ndimage.map_coordinates(wf.w, [[x1], [y1], [z1]], order=3, prefilter=false)
        # v_x = x_wind[0] * rel_turb + v_wind_height
        # v_y = y_wind[0] * rel_turb
        # v_z = z_wind[0] * rel_turb  
    else
        v_x = wf.u[x1, y1, z1] * rel_turb + v_wind_height
        v_y = wf.v[x1, y1, z1] * rel_turb
        v_z = wf.w[x1, y1, z1] * rel_turb
        return v_x, v_y, v_z
    end
    return nothing
end

#     def getWind(self, x, y, z, t, interpolate=True, rel_turb = 0.351):
#         """ Return the wind vector for a given position and time. Linear interpolation in x, y and z.
#         3.34 sec for order = 2; 37µs for order = 1
#         """
#         assert(z >= 5.)
#         if z < 10.:
#             z = 10.0
#         assert(t >= 0.)
#         while x < 0.0:
#             x += self.x_range
#         while y > self.y_max:
#             y -= self.y_range
#         while y < self.y_min:
#             y += self.y_range
#         y1 = ((y + self.y_range / 2) / GRID_STEP)
#         v_wind_height = calcWindHeight(V_WIND_GND, z)
#         # print "--->", v_wind_height
#         # sys.exit()
#         x1 = (x + t * v_wind_height) / GRID_STEP
#         while x1 > self.u.shape[0] - 1:
#             x1 -= self.u.shape[0] - 1
#         z1 = z / HEIGHT_STEP
#         if z1 > self.u.shape[2] - 1:
#             z1 = self.u.shape[2] - 1
#         if z1 < 0:
#             z1 = 0
#         if interpolate:
#             x_wind = ndimage.map_coordinates(self.u_pre, [[x1], [y1], [z1]], order=3, prefilter=False)
#             y_wind = ndimage.map_coordinates(self.v_pre, [[x1], [y1], [z1]], order=3, prefilter=False)
#             z_wind = ndimage.map_coordinates(self.w_pre, [[x1], [y1], [z1]], order=3, prefilter=False)
#             v_x, v_y, v_z = x_wind[0]*rel_turb + v_wind_height, y_wind[0]*rel_turb, z_wind[0]*rel_turb
#             # v_x, v_y, v_z = v_wind_height, 0.0, 0.0
#             v_wind = sqrt(v_x*v_x + v_y*v_y + v_z*v_z)
#             if v_wind < 1.0:
#                 print "x1, y1, z1", x1, y1, z1
#             # print " v_x, v_y, v_z",  v_x, v_y, v_z
#             return v_x, v_y, v_z
#         else:
#             vx, vy, vz = self.u[x1, y1, z1] + v_wind_height, self.v[x1, y1, z1], self.w[x1, y1, z1]
#             return vx, vy, vz


"""
    new_windfield(v_wind_gnd)

Create a new wind field object using the given ground wind velocity vector `v_wind_gnd`.

# Arguments
- `v_wind_gnd`: A scalar representing the wind velocity at ground level.

# Returns
nothing
"""
function new_windfield(am::AtmosphericModel, v_wind_gnd; prn=true)
    Random.seed!(1234) 
    prn && @info "Creating wind field for $v_wind_gnd m/s. This might take 30s or more..."
    y, x, z = create_grid(100, 4050, 500, 70)
    sigma1 = set.use_turbulence * calc_sigma1(am, v_wind_gnd)
    prn && @info "Creating wind field with sigma1 = $sigma1"
    u, v, w = create_windfield(x, y, z, sigma1=sigma1)
    # addWindSpeed(z, u)
    param = [am.set.alpha, v_wind_gnd]
    save(am, x, y, z, u, v, w, param; basename="windfield_4050_500", v_wind_gnd)
    prn && @info "Finished creating and saving wind field!"
    nothing
end

"""
    new_windfields(am::AtmosphericModel)

Create and initialize new wind fields for all ground wind speeds, defined in V_WIND_GNDS 
for the given `AtmosphericModel` instance `am`.

# Arguments
- `am::AtmosphericModel`: The atmospheric model for which wind fields are to be generated.

# Returns
- nothing

"""
function new_windfields(am::AtmosphericModel; prn=true)
    for v_wind_gnd in am.set.v_wind_gnds
        new_windfield(am, v_wind_gnd; prn)
    end
    @info "All wind fields created and saved successfully!"
    nothing
end

# if __name__ == "__main__":
#     SAVE = False # True: calculate and save new wind field; False: use saved wind field
#     if not SAVE:
#         if WIND_FIELD.valid:
#             x, y, z = WIND_FIELD.x, WIND_FIELD.y, WIND_FIELD.z
#             u, v, w = WIND_FIELD.u, WIND_FIELD.v, WIND_FIELD.w
#         else:
#             SAVE = True
#     v_wind = calcWindHeight(V_WIND_GND, 200)
#     if SAVE:
#         new_windfield(V_WIND_GND)
#         WIND_FIELD = WindField()

#     if False:
#         showGrid(x, y, z)
#     if False:
#         plotTurbulenceVsHeight(x, y, z, u, v, w)
#     if False:
#         plotWindVsY(x, y, z, u, v, w)
#     if False:
#         u2 =  u[:,1,:]
#         w2 =  w[:,1,:]
#         show2Dfield(x[:,1,:], z[:,1,:], u2, w2, scale=8.0)
#         # showGrid(X, Y, Z)
#     if True:
#         v_x, v_y, v_z = WIND_FIELD.getWind(0, 0, 197, 0)
#         print v_x, v_y, v_z
#         plotWindVsTime(0., 0., 197.)
#     if True:
#         print "sigma1", REL_TURB[2] * calcSigma1(V_WIND_GNDS[TEST])
#     if True:
#         print "V_WIND_GND at 6 m", V_WIND_GNDS[TEST]
#         print "wind at 197m:", calcWindHeight(V_WIND_GNDS[TEST], 197.0)

