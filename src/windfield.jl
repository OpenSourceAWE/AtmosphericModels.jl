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

function Base.getproperty(wf::WindField, sym::Symbol)
    if sym == :y_range
        getproperty(wf, :y_max) - getproperty(wf, :y_min) 
    elseif sym == :x_range
        getproperty(wf, :x_max) - getproperty(wf, :x_min)
    else
        getfield(wf, sym)
    end
end
"""
    WindField(am, speed; prn=true)

Load (or generate) the wind field for the ground wind `speed` and wrap it in a [`WindField`](@ref).

Throws if the field cannot be built: the settings are checked by
[`check_windfield_settings`](@ref) first, and any later failure propagates. Returning `nothing`
here, as versions before v0.3.8 did, moved the failure to the `wf !== nothing` assertion in
[`get_wind`](@ref), a stack trace away from its cause.
"""
function WindField(am, speed; prn=true)
    check_windfield_settings(am.set)
    last_speed = 0.0
    prn && @info "Loading wind field... $speed m/s"
    try
        u, v, w, param, v_wind_gnd = load_windfield(am, speed)
        x, y, z = grid_axes(am)
        valid = true
        x_max = maximum(x)
        x_min = minimum(x)
        y_max = maximum(y)
        y_min = minimum(y)
        z_max = maximum(z)
        z_min = minimum(z)
        return WindField(x_max, x_min, y_max, y_min, z_max, z_min, last_speed, valid, x, y, z, u, v, w,
                         param, v_wind_gnd)
    catch
        @error "Failed to build the wind field for $speed m/s in $(windfield_path())"
        rethrow()
    end
end

"""
    check_windfield_settings(set::Settings)

Throw an `ArgumentError` naming the offending key if `set` cannot describe a wind field.

Runs before a field is loaded or generated, so a missing or inconsistent setting is reported
against the settings file rather than surfacing later and elsewhere. `set.grid` defaulting to
`Int64[]` when the YAML does not declare it is the case that motivated this.
"""
function check_windfield_settings(set::Settings)
    length(set.grid) == 4 ||
        throw(ArgumentError("environment.grid must be [nx, ny, nz, z_min], got $(set.grid)"))
    all(>(0), set.grid[1:3]) ||
        throw(ArgumentError("environment.grid: nx, ny and nz must be positive, got $(set.grid)"))
    set.grid_step > 0 ||
        throw(ArgumentError("environment.grid_step must be positive, got $(set.grid_step)"))
    set.height_step > 0 ||
        throw(ArgumentError("environment.height_step must be positive, got $(set.height_step)"))
    for (i, axis) in enumerate(("nx", "ny"))
        set.grid[i] % set.grid_step == 0 ||
            throw(ArgumentError("environment.grid $axis = $(set.grid[i]) must be a multiple of " *
                                "grid_step = $(set.grid_step)"))
    end
    set.grid[3] % set.height_step == 0 ||
        throw(ArgumentError("environment.grid nz = $(set.grid[3]) must be a multiple of " *
                            "height_step = $(set.height_step)"))
    isempty(set.v_wind_gnds) &&
        throw(ArgumentError("environment.v_wind_gnds is empty; it lists the ground wind speeds a " *
                            "wind field can be generated for"))
    length(set.rel_turbs) == length(set.v_wind_gnds) ||
        throw(ArgumentError("environment.rel_turbs has $(length(set.rel_turbs)) entries, " *
                            "environment.v_wind_gnds has $(length(set.v_wind_gnds))"))
    for key in (:i_ref, :avg_height, :h_ref)
        getfield(set, key) > 0 ||
            throw(ArgumentError("environment.$key must be positive, got $(getfield(set, key)); " *
                                "it determines sigma1"))
    end
    nothing
end

function pfq(z)
    _₂F₁(1. /3 , 17. /6, 4. /3, z)
end

function calc_sigma1(am, v_wind_gnd)
    v_height = v_wind_gnd * calc_wind_factor(am, am.set.avg_height, Val{Int(EXP)}) 
    am.set.i_ref * (0.75 * v_height + 5.6)
end

"""
    rel_turbo(am::AtmosphericModel, v_wind = am.set.v_wind)

Find the closest relative turbulence value for a given ground wind speed.

# Arguments
- `am::AtmosphericModel`: The atmospheric model instance containing relevant parameters.
- `v_wind`: (Optional) The wind velocity to use for the calculation. Defaults to `am.set.v_wind`.

# Returns
- The computed relative turbulence value.
"""
function rel_turbo(am::AtmosphericModel, v_wind = am.set.v_wind)
    # Find the closest relative turbulence value for a given ground wind speed
    _, idx = findmin(abs.(am.set.v_wind_gnds .- v_wind))
    return am.set.rel_turbs[idx]
end

"""
    nextpow2(i)

Find 2^n that is equal to or greater than i.
"""
function nextpow2(i)
    n = 1
    while n < i
        n *= 2
    end
    n
end

const WINDFIELD_PATH = Ref("")

"""
    windfield_path()

Directory the generated `.npz` wind fields are read from and written to.

By default a `Scratch.jl` scratchspace, created on first use: the files are derived artifacts of
~1.2 GB apiece, so they belong in a cache shared by every package using this one, not in the
`data/` directory of each of them. Redirect it with [`set_windfield_path!`](@ref).
"""
function windfield_path()
    if isempty(WINDFIELD_PATH[])
        WINDFIELD_PATH[] = @get_scratch!("windfields")
    end
    WINDFIELD_PATH[]
end

"""
    set_windfield_path!(path)

Store the wind fields in `path` instead of the scratchspace, e.g. on a disk with room for them.
Pass `""` to go back to the default. The directory is created if it does not exist.
"""
function set_windfield_path!(path::AbstractString)
    !isempty(path) && mkpath(path)
    WINDFIELD_PATH[] = path
    path
end

function calc_full_name(v_wind_gnd; basename)
    joinpath(windfield_path(), basename * "_" * @sprintf("%.1f", v_wind_gnd))
end

"""
    find_windfield(set::Settings, v_wind_gnd)

Path of the stored wind field for `v_wind_gnd`, without the `.npz` suffix, or `nothing`.

Prefers the current name — the one carrying the [`param_digest`](@ref), which proves the file
matches `set` — over the two older ones, and [`windfield_path`](@ref) over `get_data_path()`, where
versions before v0.3.8 kept the files. A file found under an older name cannot be checked against
`set`; [`load`](@ref) says so.
"""
function find_windfield(set::Settings, v_wind_gnd)
    speed = @sprintf("%.1f", v_wind_gnd)
    grid_name = grid_basename(set)
    names = (calc_basename(set) * "_" * speed,
             grid_name * "_" * speed,        # before the parameter digest
             grid_name * "_1.0_" * speed)    # before use_turbulence moved to the lookup
    for name in names, dir in (windfield_path(), get_data_path())
        path = joinpath(dir, name)
        isfile(path * ".npz") && return path
    end
    return nothing
end

function save(am, u, v, w, param; basename=calc_basename(am.set), v_wind_gnd)
    fullname = calc_full_name(v_wind_gnd; basename)
    @info "Saving wind field to: $fullname.npz"
    # The x/y/z meshgrids are not stored: they are reconstructed by grid_axes from the settings.
    NPZ.npzwrite(fullname * ".npz", Dict(
        "u" => u,
        "v" => v,
        "w" => w,
        "param" => param
    ))
end

"""
    load(am::AtmosphericModel; v_wind_gnd=8.0)

Read the stored wind field for the ground wind speed `v_wind_gnd` and return `(u, v, w, param)`.

The file is located with [`find_windfield`](@ref), which also accepts the two names used before
v0.3.8. Which one was found is logged, because a name without the [`param_digest`](@ref) cannot be
checked against `am.set`. If there is no file at all, [`new_windfield`](@ref) generates one under
the current name, which takes ~30 s or more.

The `x`, `y` and `z` meshgrids are not read back even when an old file still stores them (622 MB of
1.24 GB for the default grid); the axes are rebuilt from the settings by [`grid_axes`](@ref).

Callers normally want [`load_windfield`](@ref), which picks `v_wind_gnd` from `am.set.v_wind_gnds`
and reports back which entry it used.
"""
function load(am::AtmosphericModel; v_wind_gnd=8.0)
    current = calc_full_name(v_wind_gnd; basename=calc_basename(am.set))
    fullname = find_windfield(am.set, v_wind_gnd)
    if isnothing(fullname)
        fullname = current
        @warn "Wind field file not found: $fullname.npz"
        new_windfield(am::AtmosphericModel, v_wind_gnd; prn=true)
    elseif basename(fullname) != basename(current)
        @info "Using $fullname.npz, whose name predates the parameter digest: it is assumed to " *
              "match the current settings. Delete it to regenerate one that is checked."
    elseif dirname(fullname) != windfield_path()
        @info "Using $fullname.npz; move it to $(windfield_path()) to share it between projects."
    end
    # Named vars, so a file written before v0.3.8 does not read back its 622 MB of coordinates.
    npzfile = NPZ.npzread(fullname * ".npz", ["u", "v", "w", "param"])
    return npzfile["u"], npzfile["v"], npzfile["w"], npzfile["param"]
end

"""
    load_windfield(am::AtmosphericModel, speed)

Load the wind field generated for the `am.set.v_wind_gnds` entry closest to `speed`.

Returns `(u, v, w, param, v_wind_gnd)`; the trailing `v_wind_gnd` is the grid speed that was
actually chosen, which is what [`get_wind`](@ref) needs to pick the matching `rel_turbs`.
"""
function load_windfield(am::AtmosphericModel, speed)
    # Find the index of the closest wind speed
    idx = findmin(abs.(am.set.v_wind_gnds .- speed))[2]
    v_wind_gnd = am.set.v_wind_gnds[idx]
    return (load(am; v_wind_gnd)..., v_wind_gnd)
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
    create_grid(am::AtmosphericModel)

Creates a 3D grid for the wind field model.

# Parameters
- `am`:       An instance of `AtmosphericModel` containing the settings.

# Returns Y, X and Z
Three arrays representing the generated 3D grid.
"""
function create_grid(am::AtmosphericModel)
    x_range, y_range, z_range = grid_axes(am)

    # Create meshgrid (Julia's meshgrid returns in order (x, y, z))
    X, Y, Z = ndgrid(x_range, y_range, z_range)

    return Y, X, Z  # To match the Python (y, x, z) order
end

"""
    grid_axes(am::AtmosphericModel)

The `(x, y, z)` coordinate axes of the wind field grid as ranges [m], from `am.set.grid`,
`am.set.grid_step` and `am.set.height_step`.

`x` runs downwind from zero, `y` is centered on zero and `z` starts at `am.set.grid[4]`. The
stored `.npz` holds only `u`, `v`, `w`, so these axes are rebuilt here instead of being read back.
"""
function grid_axes(am::AtmosphericModel)
    res = am.set.grid_step
    nx, ny, nz, z_min = am.set.grid
    x_range = range(0, nx, length=Int(nx/res)+1)
    y_range = range(-ny/2, ny/2, length=Int(ny/res)+1)
    z_range = range(z_min, z_min+nz, length=Int(nz/am.set.height_step)+1)
    return x_range, y_range, z_range
end

function meshgrid(x, y, z)
    X = [i for i in x, _ in y, _ in z]
    Y = [j for _ in x, j in y, _ in z]
    Z = [k for _ in x, _ in y, k in z]
    return (X, Y, Z)
end

function create_windfield(x, y, z; sigma1=nothing, gamma=3.9, ae=0.1, length_scale=33.6, rng=Random.default_rng())
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
    n_real = randn(rng, 3, 1, nx, ny, nz)
    n_imag = randn(rng, 3, 1, nx, ny, nz)
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
    get_wind(am::AtmosphericModel, x, y, z, t; upwind_dir=-π/4, interpolate=false)

Returns the wind vector at the specified position (`x`, `y`, `z`) and time `t` using the given
`AtmosphericModel` (`am`).

Uses Taylor's frozen-turbulence hypothesis: the field is advected along the mean wind direction.
The position is first rotated into the wind-aligned frame so that:
- the **along-wind** component (+ time advection) maps to whichever of the field's first two
  dimensions is **long** (the larger of `am.set.grid[1]`, `am.set.grid[2]`), avoiding short-period
  repetition during long simulations.
- the **cross-wind** component maps to the **short** dimension, so the kite stays within the
  spatial range of the field.

The long/short axis is detected from the actual array size at each call, since which of `x`/`y`
is longer depends on `am.set.grid` and differs between configurations (e.g. `[4050, 100, ...]`
vs. the `[100, 4050, ...]` default).

The stored field is scaled at lookup by `am.set.use_turbulence * rel_turbo(am, wf.v_wind_gnd)`, so
one file per ground wind speed serves every turbulence intensity and changing `use_turbulence` needs
no regeneration. The `rel_turbs` correction is taken for the speed the loaded field was generated
for, not for `am.set.v_wind`, so a field loaded at any other speed keeps its own intensity.

# Arguments
- `am::AtmosphericModel`: The atmospheric model providing environmental parameters.
- `x`, `y`, `z`: Position in the simulation (ENU) frame where the wind is evaluated. [m]
- `t`: Current simulation time. [s]
- `upwind_dir` (optional, default = `-π/4`): Direction the wind is coming FROM [rad].
  Zero is north, clockwise positive (same convention as in `calc_turbulent_wind`).
- `interpolate` (optional, default = `false`): If `true`, interpolate the turbulence trilinearly
  between the eight surrounding grid points; otherwise, use nearest-grid-point values. Interpolation
  removes the steps a kite flying through the field sees, at about 1.8x the cost of the lookup
  itself (29 ns → 51 ns per position for the vector method).

# Returns
- A tuple `(v_x, v_y, v_z)` representing the wind velocity in the wind-aligned frame [m/s],
  where `v_x` is the along-wind component (includes mean wind), `v_y` is cross-wind, `v_z` is vertical.
"""
function get_wind(am::AtmosphericModel, x, y, z, t; upwind_dir=-π/4, interpolate=false)
    @assert z >= 5.0 "Height must be at least 5 m"
    wf = am.wf
    @assert wf !== nothing "No wind field: AtmosphericModel(set) only loads one when set.use_turbulence > 0"
    wind_at(am, wf, x, y, z, t, wind_context(am, wf, upwind_dir)...; interpolate)
end

"""
    wind_context(am::AtmosphericModel, wf::WindField, upwind_dir)

Compute everything [`wind_at`](@ref) needs that does not depend on the position: the sine/cosine of
the wind direction, the turbulence scaling and the settings read out of `am.set`.

Hoisted out of the lookup so that the vector method of [`get_wind`](@ref) pays for it once per call
instead of once per position — `rel_turbo` in particular allocates.
"""
@inline function wind_context(am::AtmosphericModel, wf::WindField, upwind_dir)
    rel_turb = am.set.use_turbulence * rel_turbo(am, wf.v_wind_gnd)
    # wind_dir: direction the wind is blowing TO, measured from +x (East) axis, CCW.
    wind_dir = -upwind_dir - pi/2
    return (cos(wind_dir), sin(wind_dir), rel_turb, Int64(am.set.profile_law), am.set.v_wind,
            am.set.grid_step, am.set.height_step)
end

"""
    interp_bounds_periodic(idx, n)

Split the fractional grid index `idx`, already reduced into `[0, n-1]`, into the pair of neighbouring
one-based indices and the weight of the upper one. The upper neighbour wraps around with the same
period of `n - 1` grid steps that the nearest-grid-point lookup uses.
"""
@inline function interp_bounds_periodic(idx, n)
    lower = floor(idx)
    i_lo = Int(lower) + 1
    i_hi = i_lo + 1
    if i_hi > n
        i_hi -= n - 1
    end
    return i_lo, i_hi, idx - lower
end

"""
    interp_bounds_clamped(idx, n)

As [`interp_bounds_periodic`](@ref), but for the vertical axis, which is not periodic: at the top of
the grid the upper neighbour is clamped to the last layer.
"""
@inline function interp_bounds_clamped(idx, n)
    lower = floor(idx)
    i_lo = Int(lower) + 1
    return i_lo, min(i_lo + 1, n), idx - lower
end

"""
    trilinear(a, i_lo, i_hi, j_lo, j_hi, k_lo, k_hi, fi, fj, fk)

Trilinear interpolation of the 3D array `a` between the eight grid points given by the index pairs,
with `fi`, `fj`, `fk` the weights of the upper index in each dimension.
"""
@inline function trilinear(a, i_lo, i_hi, j_lo, j_hi, k_lo, k_hi, fi, fj, fk)
    a00 = a[i_lo, j_lo, k_lo] * (1 - fi) + a[i_hi, j_lo, k_lo] * fi
    a10 = a[i_lo, j_hi, k_lo] * (1 - fi) + a[i_hi, j_hi, k_lo] * fi
    a01 = a[i_lo, j_lo, k_hi] * (1 - fi) + a[i_hi, j_lo, k_hi] * fi
    a11 = a[i_lo, j_hi, k_hi] * (1 - fi) + a[i_hi, j_hi, k_hi] * fi
    a0 = a00 * (1 - fj) + a10 * fj
    a1 = a01 * (1 - fj) + a11 * fj
    return a0 * (1 - fk) + a1 * fk
end

"""
    wind_at(am, wf, x, y, z, t, cos_dir, sin_dir, rel_turb, profile_law, v_wind, grid_step,
            height_step; interpolate=false)

Look up the wind vector at one position; the trailing arguments come from [`wind_context`](@ref).

With `interpolate=false` the turbulence is read at the nearest grid point, with `interpolate=true` it
is interpolated trilinearly between the eight surrounding ones. The mean wind is a smooth function of
the height either way, so only the turbulence is affected.

Returns the tuple `(v_x, v_y, v_z)` in the wind-aligned frame, see [`get_wind`](@ref).
"""
@inline function wind_at(am::AtmosphericModel, wf::WindField, x, y, z, t, cos_dir, sin_dir, rel_turb,
                         profile_law, v_wind, grid_step, height_step; interpolate=false)
    @assert z >= 5.0 "Height must be at least 5 m"
    @assert t >= 0.0 "Time must be non-negative"
    if z < 10.0
        z = 10.0
    end

    # Rotate (x, y) from simulation (ENU) frame into wind-aligned frame.
    along = x * cos_dir + y * sin_dir   # along-wind component
    cross = -x * sin_dir + y * cos_dir  # cross-wind component

    v_wind_height = v_wind * calc_wind_factor(am, z, profile_law)

    n1 = size(wf.u, 1)
    n2 = size(wf.u, 2)
    nz = size(wf.u, 3)
    dim1_is_long = n1 >= n2
    nlong = dim1_is_long ? n1 : n2
    nshort = dim1_is_long ? n2 : n1

    # Along-wind + Taylor advection → long field dimension (avoids short-period repetition)
    along_idx = (along + t * v_wind_height) / grid_step
    while along_idx > nlong - 1
        along_idx -= nlong - 1
    end
    while along_idx < 0
        along_idx += nlong - 1
    end

    # Cross-wind → short field dimension (kite stays within spatial range)
    cross_idx = cross / grid_step
    while cross_idx > nshort - 1
        cross_idx -= nshort - 1
    end
    while cross_idx < 0
        cross_idx += nshort - 1
    end

    z1 = z / height_step
    if z1 > nz - 1
        z1 = nz - 1
    elseif z1 < 0
        z1 = 0
    end

    if interpolate
        long_lo, long_hi, f_long = interp_bounds_periodic(along_idx, nlong)
        short_lo, short_hi, f_short = interp_bounds_periodic(cross_idx, nshort)
        k_lo, k_hi, f_z = interp_bounds_clamped(z1, nz)
        i_lo, i_hi, f_i = dim1_is_long ? (long_lo, long_hi, f_long) : (short_lo, short_hi, f_short)
        j_lo, j_hi, f_j = dim1_is_long ? (short_lo, short_hi, f_short) : (long_lo, long_hi, f_long)
        v_x = trilinear(wf.u, i_lo, i_hi, j_lo, j_hi, k_lo, k_hi, f_i, f_j, f_z) * rel_turb + v_wind_height
        v_y = trilinear(wf.v, i_lo, i_hi, j_lo, j_hi, k_lo, k_hi, f_i, f_j, f_z) * rel_turb
        v_z = trilinear(wf.w, i_lo, i_hi, j_lo, j_hi, k_lo, k_hi, f_i, f_j, f_z) * rel_turb
        return v_x, v_y, v_z
    end

    i_long = Int(round(along_idx)) + 1
    i_short = Int(round(cross_idx)) + 1
    k = Int(round(z1)) + 1

    i = dim1_is_long ? i_long : i_short
    j = dim1_is_long ? i_short : i_long

    v_x = wf.u[i, j, k] * rel_turb + v_wind_height
    v_y = wf.v[i, j, k] * rel_turb
    v_z = wf.w[i, j, k] * rel_turb
    return v_x, v_y, v_z
end

"""
    get_wind(am::AtmosphericModel, positions::AbstractVector, t; upwind_dir=-π/4, interpolate=false)

Return the wind vectors at all `positions` at time `t` as a `Vector{SVec3}`.

Faster than calling the scalar [`get_wind`](@ref) in a loop: the turbulence scaling
(`rel_turbo`, which allocates), the sine/cosine of the wind direction and the settings lookups are
computed once per call instead of once per position.

# Arguments
- `am::AtmosphericModel`: The atmospheric model providing environmental parameters.
- `positions`: Vector of 3D positions in the simulation (ENU) frame, e.g. `Vector{SVec3}`; any
  vector of indexable, 3-element positions works. The height `positions[i][3]` must be >= 5 m.
- `t`: Current simulation time. [s]
- `upwind_dir` (optional, default = `-π/4`): Direction the wind is coming FROM [rad].
- `interpolate` (optional, default = `false`): Interpolate between grid points, see the scalar
  method of [`get_wind`](@ref).

# Returns
- A `Vector{SVec3}` of wind velocities in the wind-aligned frame [m/s], one per position; the
  components are `(v_x, v_y, v_z)` as returned by the scalar method.

See also [`get_wind!`](@ref), which writes into a pre-allocated result vector.
"""
function get_wind(am::AtmosphericModel, positions::AbstractVector, t; upwind_dir=-π/4, interpolate=false)
    res = Vector{SVec3}(undef, length(positions))
    get_wind!(res, am, positions, t; upwind_dir, interpolate)
end

"""
    get_wind!(res::AbstractVector{SVec3}, am::AtmosphericModel, positions::AbstractVector, t;
              upwind_dir=-π/4, interpolate=false)

In-place version of the vector method of [`get_wind`](@ref): write the wind vectors at `positions`
into `res` and return `res`. `res` must have the same length as `positions`.
"""
function get_wind!(res::AbstractVector{SVec3}, am::AtmosphericModel, positions::AbstractVector, t;
                   upwind_dir=-π/4, interpolate=false)
    length(res) == length(positions) ||
        throw(DimensionMismatch("res has $(length(res)) entries, positions $(length(positions))"))
    wf = am.wf
    @assert wf !== nothing "No wind field: AtmosphericModel(set) only loads one when set.use_turbulence > 0"
    ctx = wind_context(am, wf, upwind_dir)
    # No @inbounds: the field lookup in wind_at keeps the same bounds checks as the scalar method.
    for i in eachindex(positions, res)
        pos = positions[i]
        res[i] = SVec3(wind_at(am, wf, pos[1], pos[2], pos[3], t, ctx...; interpolate))
    end
    return res
end

"""
    calc_turbulent_wind(am::AtmosphericModel, pos, t; upwind_dir=-π/4, interpolate=false)

Calculate the wind velocity vectors at the kite and at the mid-tether point, in the ENU
simulation frame.

When `am.set.use_turbulence == 0`, the mean wind for the configured `am.set.profile_law` is
returned. Otherwise the turbulent wind vectors are looked up from the pre-computed wind field
via [`get_wind`](@ref) and rotated from the wind-aligned frame into the simulation frame.

# Arguments
- `am::AtmosphericModel`: atmospheric model; the settings are read from `am.set`.
- `pos`: 3D position of the kite [m]; `pos[3]` is the height, clamped to $(MIN_KITE_HEIGHT) m minimum.
- `t`: current simulation time [s].
- `upwind_dir` (optional, default = `-π/4`): direction the wind is coming FROM [rad].
  Zero is north, clockwise positive (same convention as in [`get_wind`](@ref)).
- `interpolate` (optional, default = `false`): interpolate between grid points, see [`get_wind`](@ref).

# Returns
A tuple `(v_wind, v_wind_tether)` of `SVec3` in the ENU frame [m/s]:
- `v_wind`: wind velocity at the kite position.
- `v_wind_tether`: wind velocity at half the kite position `(0.5x, 0.5y, 0.5z)`, with the height
  clamped to $(MIN_TETHER_HEIGHT) m minimum.
"""
function calc_turbulent_wind(am::AtmosphericModel, pos, t; upwind_dir=-π/4, interpolate=false)
    wind_dir = -upwind_dir - pi/2
    height = max(pos[3], MIN_KITE_HEIGHT)
    rotate_wind(wx, wy, wz) = SVec3(wx * cos(wind_dir) - wy * sin(wind_dir),
                                    wx * sin(wind_dir) + wy * cos(wind_dir),
                                    wz)
    if am.set.use_turbulence == 0
        v_wind_gnd = am.set.v_wind
        mean_wind(h) = SVec3(v_wind_gnd * calc_wind_factor(am, h) * cos(wind_dir),
                             v_wind_gnd * calc_wind_factor(am, h) * sin(wind_dir),
                             0.0)
        return mean_wind(height), mean_wind(height / 2.0)
    end
    v_wind = rotate_wind(get_wind(am, pos[1], pos[2], height, t; upwind_dir, interpolate)...)
    tether_height = max(0.5 * height, MIN_TETHER_HEIGHT)
    v_wind_tether = rotate_wind(get_wind(am, 0.5 * pos[1], 0.5 * pos[2], tether_height, t;
                                         upwind_dir, interpolate)...)
    return v_wind, v_wind_tether
end

"""
    new_windfield(am::AtmosphericModel, v_wind_gnd; prn=true)

Create a new wind field file using the given, scalar ground wind velocity `v_wind_gnd`.

The field is stored at the reference intensity (`sigma1 = calc_sigma1(am, v_wind_gnd)`);
`am.set.use_turbulence` and `rel_turbo(am)` are applied when it is read by [`get_wind`](@ref).

# Parameters
- `am::AtmosphericModel`: The atmospheric model for which the wind field is created.
- `v_wind_gnd`: A scalar representing the wind velocity at ground level.
- `prn`: Optional boolean flag to control printing of progress messages (default is `true`).

# Returns
- nothing
"""
function new_windfield(am::AtmosphericModel, v_wind_gnd; prn=true)
    check_windfield_settings(am.set)
    prn && @info "Creating wind field for $v_wind_gnd m/s. This might take 30s or more..."
    y, x, z = create_grid(am)
    sigma1 = calc_sigma1(am, v_wind_gnd)
    u, v, w = create_windfield(x, y, z, sigma1=sigma1, rng=StableRNG(1234))
    param = [am.set.alpha, v_wind_gnd]
    save(am, u, v, w, param; basename=calc_basename(am.set), v_wind_gnd)
    prn && @info "Finished creating and saving wind field!"
    nothing
end

"""
    grid_basename(set::Settings)

File name prefix of a wind field, from `set.grid` alone. The name a version before v0.3.8 gave the
file; today [`calc_basename`](@ref) appends the [`param_digest`](@ref) to it.
"""
function grid_basename(set::Settings)
    "windfield_$(string(set.grid[1]))_$(string(set.grid[2]))_$(string(set.grid[3]))_$(string(set.grid[4]))"
end

"""
    param_digest(set::Settings)

Eight hex digits of a SHA-256 over the settings that change the generated field but are not in its
file name: `grid_step` and `height_step`, which resolve the grid, and `i_ref`, `alpha`, `avg_height`
and `h_ref`, which enter `calc_sigma1`. `grid` and the ground wind speed are in the name already.

`profile_law` and `z0` are deliberately not in it: neither reaches the generator, since
`calc_sigma1` evaluates the wind profile as `EXP` whatever the setting says. `use_turbulence` is
not either — it scales the field at lookup, see [`get_wind`](@ref).
"""
function param_digest(set::Settings)
    key = string("grid_step=", set.grid_step, ";height_step=", set.height_step,
                 ";i_ref=", set.i_ref, ";alpha=", set.alpha,
                 ";avg_height=", set.avg_height, ";h_ref=", set.h_ref)
    bytes2hex(sha256(key))[1:8]
end

"""
    calc_basename(set::Settings)

File name prefix of a wind field: [`grid_basename`](@ref) plus the [`param_digest`](@ref), so that
changing any setting the field depends on names a different file instead of silently loading the
old one.
"""
calc_basename(set::Settings) = grid_basename(set) * "_" * param_digest(set)

"""
    new_windfields(am::AtmosphericModel; prn=true)

Create and initialize new wind fields for all ground wind speeds, defined in `am.set.v_wind_gnds` and save them
for the given `AtmosphericModel` instance `am`.

# Arguments
- `am::AtmosphericModel`: The atmospheric model for which wind fields are to be generated.
- `prn`: Optional boolean flag to control printing of progress messages (default is `true`).

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
