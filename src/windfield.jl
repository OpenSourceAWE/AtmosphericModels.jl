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
const REL_SIGMA  =  1.0    # turbulence relative to the IEC model
const V_WIND_GND  = 8.0    # Default value, change as needed
const GRID_STEP   = 2.0    # Resolution of the grid in x and y direction in meters
const HEIGHT_STEP = 2.0    # Resolution in z direction in meters

function pfq(z)
    _₂F₁(1. /3 , 17. /6, 4. /3, z)
end

function calc_sigma1(am, v_wind_gnd)
    v_height = v_wind_gnd * calc_wind_factor(am, am.set.avg_height, Val{Int(EXP)})
    am.set.i_ref * (0.75 * v_height + 5.6)
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

function calcFullName(v_wind_gnd; basename="windfield_4050_500", rel_sigma=1.0)
    path = get_data_path() * "/"
    name = basename * "_" * @sprintf("%.1f", rel_sigma)
    name *= "_" * @sprintf("%.1f", v_wind_gnd)
    return path * name
end

function save(x, y, z, u, v, w, param; basename="windfield_4050_500", v_wind_gnd=V_WIND_GND)
    fullname = calcFullName(v_wind_gnd; basename)
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
    fullname = calcFullName(v_wind_gnd, basename=basename)
    npzfile = NPZ.npzread(fullname * ".npz")
    return npzfile["x"], npzfile["y"], npzfile["z"], npzfile["u"], npzfile["v"], npzfile["w"], npzfile["param"]
end

V_WIND_GNDS = [3.483, 5.324, 8.163]

function load_windfield(speed)
    # Find the index of the closest wind speed
    idx = findmin(abs.(V_WIND_GNDS .- speed))[2]
    return load(v_wind_gnd = V_WIND_GNDS[idx])
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

# def show2Dfield(X, Y, U, V, scale=1.0, width=0.07):
#     """
#     x, y: position of the wind velocity vectors
#     u, v: wind velocity vector components
#     """
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#     ax.quiver(X, Y, U, V, angles='uv', units='xy', scale=scale, width=width)

# def plotTurbulenceVsHeight(X, Y, Z, U, V, W):
#     fig = plt.figure()
#     height_min = np.min(Z)
#     height_max = np.max(Z)
#     print "height range: ", height_min, height_max
#     Height, XX, YY, ZZ = [], [], [], []
#     for height in np.linspace(height_min, height_max, int(round((height_max-height_min) / HEIGHT_STEP))+1):
#         Height.append(height)
#         v_wind = calcWindHeight(V_WIND_GND, height)
#         height_index = int(round((height-height_min) / HEIGHT_STEP))
#         u2 =  U[:,:,height_index]; XX.append(u2.std() / v_wind * 100.0)
#         v2 =  V[:,:,height_index]; YY.append(v2.std() / v_wind * 100.0)
#         w2 =  W[:,:,height_index]; ZZ.append(w2.std() / v_wind * 100.0)
#     plt.plot(Height, XX, label = 'rel. turbulence x [%]', color="black")
#     plt.plot(Height, YY, label = 'rel. turbulence y [%]', color="blue")
#     plt.plot(Height, ZZ, label = 'rel. turbulence z [%]', color="green")
#     plt.legend(loc='upper right')
#     plt.gca().set_ylabel('Rel. turbulence [%]')
#     plt.gca().set_xlabel('Height [m]')

function meshgrid(x, y, z)
    X = [i for i in x, _ in y, _ in z]
    Y = [j for _ in x, j in y, _ in z]
    Z = [k for _ in x, _ in y, k in z]
    return (X, Y, Z)
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
function create_windfield(x::AbstractArray, y::AbstractArray, z::AbstractArray;
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
    min_height = minimum(z)
    max_height = maximum(z)
    println("Minimal height: ", min_height)
    println("Maximal height: ", max_height)
    for i in axes(z, 3)
        height = z[1, 1, i]
        v_wind = am.set.v_wind * calc_wind_factor(am, height, am.set.profile_law)
        u[:, :, i] .+= v_wind
    end
end

const SRL = StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}

Base.@kwdef struct WindField
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
function WindField(speed)
    try
        last_speed = 0.0
        println("Loading wind field... $speed m/s")
        x, y, z, u, v, w, param = load_windfield(speed)
        valid = true
        return WindField(last_speed, valid, x, y, z, u, v, w, param)
    catch
        @error "Error reading wind field!"
        return nothing
    end
end


# class WindField(object):
#     def __init__(self, speed = V_WIND_GND):
#         try:
#             self.last_speed = 0.0
#             self.load(speed)
#             self.valid = True
#         except:
#            self.valid = False
#            print "Error reading wind field!"

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

function Base.getproperty(wf::WindField, sym::Symbol)
    if sym == :x_max
        maximum(getproperty(wf, :x))
    elseif sym == :x_min
        minimum(getproperty(wf, :x))
    elseif sym == :y_max
        maximum(getproperty(wf, :y))
    elseif sym == :y_min
        minimum(getproperty(wf, :y))
    elseif sym == :y_range
        getproperty(wf, :y_max) - getproperty(wf, :y_min) 
    elseif sym == :x_range
        getproperty(wf, :x_max) - getproperty(wf, :x_min)
    else
        getfield(wf, sym)
    end
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

# WIND_FIELD = WindField()

# def plotWindVsTime(x=0.0, y=0.0, z=197.3, rel_turb = REL_TURB[TEST]):
#     print "Relative turbulence: ", rel_turb
#     fig = plt.figure()
#     TIME, v_wind_x, v_wind_norm = [], [], []
#     for t in np.linspace(0.0, 600.0, 600*20):
#         TIME.append(t)
#         v_x, v_y, v_z = WIND_FIELD.getWind(x, y, z, t, rel_turb = rel_turb)
#         # print " v_x, v_y, v_z",  v_x, v_y, v_z
#         # sys.exit()
#         v_wind = sqrt(v_x*v_x + v_y*v_y + v_z*v_z)
#         v_wind_norm.append(v_wind)
#         v_wind_x.append(v_x)
#         if v_wind < 0.1:
#             print "Error for x, y, z, t: ", x, y, z, t
#         # print t, v_x
#     v_wind_x = np.array(v_wind_x)
#     su = np.std(v_wind_x)
#     mean =  v_wind_x.mean()
#     print "Mean wind x, standart deviation, turbulence intensity [%]: ", form(mean), form(su), \
#                                                                          form(su/mean * 100.0)
#     plt.plot(TIME, v_wind_x, label = 'Abs. wind speed at 197.3 m [m/s]', color="black")
#     # plt.legend(loc='upper right')
#     plt.gca().set_xlabel('Time [s]')
#     plt.gca().set_ylabel('Abs. wind speed at 197.3 m height [m/s]')

# def plotWindVsY(X, Y, Z, U, V, W, x=200, z=200.0):
#     fig = plt.figure()
#     y_min = np.min(Y)
#     y_max = np.max(Y)
#     print "y range: ", y_min, y_max
#     YY, TURB_1, TURB_2, TURB_3 = [], [], [], []
#     if False:
#         for y in np.linspace(y_min, y_max, int(round((y_max-y_min) / GRID_STEP)) + 1):
#             YY.append(y)
#             u2_1 = U[x, y, (z-100.) / HEIGHT_STEP]; TURB_1.append(u2_1)
#             u2_2 = U[x, y, z        / HEIGHT_STEP]; TURB_2.append(u2_2)
#             u2_3 = U[x, y, (z+100.) / HEIGHT_STEP]; TURB_3.append(u2_3)
#     else:
#         t = 0.0
#         for y in np.linspace(2*y_min, 2*y_max, 1000):
#             YY.append(y)
#             u2_1, v, w = WIND_FIELD.getWind(x, y, z-100., t) #U[x, y, (z-100.) / HEIGHT_STEP];
#             TURB_1.append(u2_1)
#             u2_2, v, w= WIND_FIELD.getWind(x, y, z, t) #U[x, y, z        / HEIGHT_STEP];
#             TURB_2.append(u2_2)
#             u2_3, v, w = WIND_FIELD.getWind(x, y, z+100., t) #U[x, y, (z+100.) / HEIGHT_STEP];
#             TURB_3.append(u2_3)
#     #YY.extend((y_max - y_min + GRID_STEP) + np.array(YY))
#     #TURB_1.extend(TURB_1); TURB_2.extend(TURB_2); TURB_3.extend(TURB_3)
#     plt.plot(YY, TURB_1, label = 'abs. wind x at 100m [m/s]', color="black")
#     plt.plot(YY, TURB_2, label = 'abs. wind x at 200m [m/s]', color="blue")
#     plt.plot(YY, TURB_3, label = 'abs. wind x at 300m [m/s]', color="red")
#     plt.legend(loc='upper right')
#     plt.gca().set_xlabel('Y position [m]')

"""
    new_windfield(v_wind_gnd)

Create a new wind field object using the given ground wind velocity vector `v_wind_gnd`.

# Arguments
- `v_wind_gnd`: A scalar representing the wind velocity at ground level.

# Returns
nothing
"""
function new_windfield(am::AtmosphericModel, v_wind_gnd)
    @info "Creating wind field. This might take 30s or more..."
    y, x, z = create_grid(100, 4050, 500, 70)
    sigma1 = REL_SIGMA * calc_sigma1(am, v_wind_gnd)
    u, v, w = create_windfield(x, y, z, sigma1=sigma1)
    # addWindSpeed(z, u)
    param = [am.set.alpha, v_wind_gnd]
    save(x, y, z, u, v, w, param; basename="windfield_4050_500", v_wind_gnd=v_wind_gnd)
    @info "Finished creating and saving wind field!"
    nothing
end

# def new_windfields():
#     """
#     Create and save new wind fields for all ground wind speeds, defined in V_WIND_GNDS.
#     """

#     for v_wind_gnd in V_WIND_GNDS:
#         print "Creating wind field. This might take 30s or more..."
#         y, x, z = create_grid(100, 4050, 500, 70)
#         sigma1 = REL_SIGMA * calcSigma1(v_wind_gnd)
#         u, v, w = create_windfield(x, y, z, sigma1=sigma1)
#         param = np.array((ALPHA, v_wind_gnd))
#         save(x, y, z, u, v, w, param, v_wind_gnd=v_wind_gnd)
#         print "Finished creating and saving wind field!"
#         del y, x, z, u, v, w

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

