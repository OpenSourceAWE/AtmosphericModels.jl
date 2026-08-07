module AtmosphericModels

using KiteUtils
using KiteUtils: SVec3
using HypergeometricFunctions:_₂F₁
using NPZ, Printf
using FFTW, LinearAlgebra, Random, Statistics

export AtmosphericModel, ProfileLaw, WindField, EXP, LOG, EXPLOG, CONSTANT
export CUSTOM_LOG, CUSTOM_EXP, CUSTOM_JET
export clear, calc_rho, calc_wind_factor, rel_turbo, custom_log, custom_exp, custom_jet

export new_windfield, new_windfields, get_wind, calc_turbulent_wind

const ABS_ZERO = -273.15
const SRL = StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
const MIN_KITE_HEIGHT = 6.0
const MIN_TETHER_HEIGHT = 5.0

"""
    struct WindField

Struct that is storing a 3D model of wind vectors of the atmosphere. The Fields
x, y and z store the grid coordinates, the fields u, v and w the wind turbulence
vectors. 

# Fields
- x_max::Float64 = NaN
- x_min::Float64 = NaN
- y_max::Float64 = NaN
- y_min::Float64 = NaN
- z_max::Float64 = NaN
- z_min::Float64 = NaN
- last_speed::Float64 = 0.0
- valid::Bool = false
- x::Union{SRL, Array{Float64, 3}}
- y::Union{SRL, Array{Float64, 3}}
- z::Union{SRL, Array{Float64, 3}}
- u::Array{Float64, 3}
- v::Array{Float64, 3}
- w::Array{Float64, 3}
- param::Vector{Float64} = [0, 0] # [alpha, `v_wind_gnd`]
"""
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

"""
    mutable struct AtmosphericModel

Struct that is storing the settings and the state of the atmosphere.

# Fields
- set::Settings: The Settings struct
- `rho_zero_temp`
- wf::Union{WindField, Nothing}: The 3D [`WindField`](@ref) or `nothing`
- `jet_cache`: cached `(heights, speeds, coeffs)` of the last [`custom_jet`](@ref) fit, or `nothing`
"""
Base.@kwdef mutable struct AtmosphericModel
    set::Settings
    rho_zero_temp::Float64 = (15.0 - ABS_ZERO) / (set.temp_ref - ABS_ZERO) * set.rho_0
    wf::Union{WindField, Nothing} = nothing
    jet_cache::Union{Nothing, Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}} = nothing
end

"""
    AtmosphericModel(set::Settings; nowindfield::Bool=false)

Constructs an `AtmosphericModel` using the provided `Settings`.

# Arguments
- `set::Settings`: The settings object containing configuration parameters for the atmospheric model.
- `nowindfield::Bool=false`: Optional keyword argument. If `true`, the wind field will not be loaded.

# Returns
- An instance of `AtmosphericModel` configured according to the provided settings.
"""
function AtmosphericModel(set::Settings; nowindfield::Bool=false) 
    am = AtmosphericModel(set=set)
    if set.use_turbulence > 0 && !nowindfield
        am.wf = WindField(am, am.set.v_wind)
    end
    am
end

const AM = AtmosphericModel

"""
    clear(s::AM)

Clears or resets the state of the given `AM` (Atmospheric Model) instance `s`.

# Arguments
- `s::AM`: An instance of the `AM` (Atmospheric Model) struct containing atmospheric parameters.

# Returns
- nothing
"""
function clear(s::AM)
     s.rho_zero_temp       = (15.0 - ABS_ZERO) / (s.set.temp_ref - ABS_ZERO) * s.set.rho_0
    nothing
end

@inline function fastexp(x)
  y = 1 + x / 1024
  y *= y; y *= y; y *= y; y *= y; y *= y  
  y *= y; y *= y; y *= y; y *= y; y *= y 
  y
end

"""
    calc_rho(s::AM, height)

Calculates the air density at a given height above ground level.

# Arguments
- `s::AM`: An instance of the `AM` (Atmospheric Model) struct containing atmospheric parameters.
- `height`: The height above ground level (in meters) at which to calculate the air density.

# Returns
- The air density at the specified height (in kg/m³).

# Notes
- The calculation assumes an exponential decrease of air density with altitude.
- `s.rho_zero_temp` is the reference air density at ground level.
- `s.set.height_gnd` is the ground height offset.
- The scale height used is 8550.0 meters.
"""
calc_rho(s::AM, height) = s.rho_zero_temp * fastexp(-(height+s.set.height_gnd) / 8550.0)

"""
    @enum ProfileLaw CONSTANT=0 EXP=1 LOG=2 EXPLOG=3 CUSTOM_LOG=4 CUSTOM_EXP=5 CUSTOM_JET=6

Enumeration to describe the wind profile law that is used.
""" ProfileLaw

@enum ProfileLaw begin
    CONSTANT=0
    EXP=1
    LOG=2
    EXPLOG=3
    CUSTOM_LOG=4
    CUSTOM_EXP=5
    CUSTOM_JET=6
end

@doc """
    CONSTANT::ProfileLaw

Constant wind profile.

See also [`ProfileLaw`](@ref).
""" CONSTANT
@doc """
    EXP::ProfileLaw

Exponential wind profile.

See also [`ProfileLaw`](@ref).
""" EXP
@doc """
    LOG::ProfileLaw

Logarithmic wind profile.

See also [`ProfileLaw`](@ref).
""" LOG
@doc """
    EXPLOG::ProfileLaw

A linear combination of exponential and logarithmic wind profile to match a specific site.

See also [`ProfileLaw`](@ref).
""" EXPLOG
@doc """
    CUSTOM_LOG::ProfileLaw

Logarithmic wind profile fitted to `am.set.heights`/`am.set.speeds` by least squares, see
[`custom_log`](@ref).

See also [`ProfileLaw`](@ref).
""" CUSTOM_LOG
@doc """
    CUSTOM_EXP::ProfileLaw

Power-law wind profile fitted to `am.set.heights`/`am.set.speeds` by least squares, see
[`custom_exp`](@ref).

See also [`ProfileLaw`](@ref).
""" CUSTOM_EXP
@doc """
    CUSTOM_JET::ProfileLaw

Low-level jet wind profile, a power-law background plus a Gaussian bump, fitted to
`am.set.heights`/`am.set.speeds` by nonlinear least squares, see [`custom_jet`](@ref).

See also [`ProfileLaw`](@ref).
""" CUSTOM_JET


# Calculate the wind speed at a given height and reference height.
@inline function calc_wind_factor1(s::AM, height);  exp(s.set.alpha * log(height/s.set.h_ref)); end
@inline function calc_wind_factor(s::AM, height, ::Type{Val{0}}); 1.0; end
@inline function calc_wind_factor(s::AM, height, ::Type{Val{1}})
    calc_wind_factor1(s, height)
end 

@inline function calc_wind_factor2(s::AM, height);  log(height / s.set.z0) / log(s.set.h_ref / s.set.z0); end
@inline function calc_wind_factor(s::AM, height, ::Type{Val{2}})
    calc_wind_factor2(s, height)
end

@inline function calc_wind_factor(s::AM, height, ::Type{Val{3}}); calc_wind_factor3(s, height); end
@inline function calc_wind_factor3(s::AM, height)
    K = 1.0
    log1 = log(height / s.set.z0) / log(s.set.h_ref / s.set.z0)
    exp1 = exp(s.set.alpha * log(height/s.set.h_ref))
    log1 +  K * (log1 - exp1)
end

# Ordinary least squares fit of y = slope*x + intercept, returns (slope, intercept).
function linreg(x::AbstractVector, y::AbstractVector)
    xm = sum(x) / length(x)
    ym = sum(y) / length(y)
    sxx = sum((x .- xm) .^ 2)
    sxy = sum((x .- xm) .* (y .- ym))
    slope = sxy / sxx
    slope, ym - slope * xm
end

function check_heights_speeds(heights, speeds)
    length(heights) == length(speeds) ||
        throw(ArgumentError("heights and speeds must have the same length"))
    length(heights) >= 2 ||
        throw(ArgumentError("need at least two height/speed pairs to fit a profile"))
end

"""
    custom_log(heights, speeds, height)

Evaluate a logarithmic wind profile `u(z) = a * log(z) + b` at `height`, with `a` and `b`
fitted to `heights`/`speeds` by ordinary least squares.
"""
function custom_log(heights::AbstractVector, speeds::AbstractVector, height)
    check_heights_speeds(heights, speeds)
    slope, intercept = linreg(log.(heights), speeds)
    slope * log(height) + intercept
end

"""
    custom_exp(heights, speeds, height)

Evaluate a power-law wind profile `u(z) = c * z^a` at `height`, with `a` and `c` fitted to
`heights`/`speeds` by ordinary least squares in log-log space.
"""
function custom_exp(heights::AbstractVector, speeds::AbstractVector, height)
    check_heights_speeds(heights, speeds)
    slope, intercept = linreg(log.(heights), log.(speeds))
    exp(intercept) * height^slope
end

@inline function calc_wind_factor4(s::AM, height); custom_log(s.set.heights, s.set.speeds, height) / s.set.v_wind; end
@inline function calc_wind_factor(s::AM, height, ::Type{Val{4}}); calc_wind_factor4(s, height); end

@inline function calc_wind_factor5(s::AM, height); custom_exp(s.set.heights, s.set.speeds, height) / s.set.v_wind; end
@inline function calc_wind_factor(s::AM, height, ::Type{Val{5}}); calc_wind_factor5(s, height); end

jet_model(z, p) = p[1] * z^p[2] + p[3] * exp(-(z - p[4])^2 / (2 * p[5]^2))

function jet_jacobian_residual(heights, speeds, p)
    c, a, U_J, z_c, sigma = p
    n = length(heights)
    R = zeros(Float64, n)
    J = zeros(Float64, n, 5)
    for i in 1:n
        z = heights[i]
        za = z^a
        g = exp(-(z - z_c)^2 / (2 * sigma^2))
        R[i] = c * za + U_J * g - speeds[i]
        J[i, 1] = za
        J[i, 2] = c * za * log(z)
        J[i, 3] = g
        J[i, 4] = U_J * g * (z - z_c) / sigma^2
        J[i, 5] = U_J * g * (z - z_c)^2 / sigma^3
    end
    R, J
end

# Levenberg-Marquardt fit of jet_model's 5 coefficients [c, a, U_J, z_c, sigma] to heights/speeds.
function fit_jet(heights::AbstractVector, speeds::AbstractVector; max_iter=200, tol=1e-10)
    a0, c0 = linreg(log.(heights), log.(speeds))
    c0 = exp(c0)
    resid0 = speeds .- c0 .* heights .^ a0
    idx = argmax(abs.(resid0))
    p = Float64[c0, a0, resid0[idx], heights[idx], max((maximum(heights) - minimum(heights)) / 4, eps())]

    lambda = 1e-3
    R, J = jet_jacobian_residual(heights, speeds, p)
    cost = sum(abs2, R)
    for _ in 1:max_iter
        JtJ = J' * J
        delta = (JtJ + lambda * Diagonal(JtJ)) \ (-(J' * R))
        p_new = p + delta
        R_new, J_new = jet_jacobian_residual(heights, speeds, p_new)
        cost_new = sum(abs2, R_new)
        if cost_new < cost
            p, R, J, cost = p_new, R_new, J_new, cost_new
            lambda = max(lambda / 10, 1e-12)
            norm(delta) < tol && break
        else
            lambda *= 10
        end
    end
    p[5] = abs(p[5])
    p
end

function check_jet_data(heights, speeds)
    check_heights_speeds(heights, speeds)
    length(heights) >= 5 ||
        throw(ArgumentError("need at least five height/speed pairs to fit a jet profile"))
end

"""
    custom_jet(heights, speeds, height)

Evaluate a power-law background profile with a superimposed Gaussian jet,
`u(z) = c * z^a + U_J * exp(-(z-z_c)^2/(2*sigma^2))`, at `height`. All five coefficients
are fitted jointly to `heights`/`speeds` by nonlinear least squares (Levenberg-Marquardt).
"""
function custom_jet(heights::AbstractVector, speeds::AbstractVector, height)
    check_jet_data(heights, speeds)
    jet_model(height, fit_jet(heights, speeds))
end

# Fit CUSTOM_JET's coefficients, reusing s.jet_cache while am.set.heights/speeds are unchanged.
function jet_coeffs(s::AM)
    heights, speeds = s.set.heights, s.set.speeds
    cache = s.jet_cache
    if cache !== nothing && cache[1] == heights && cache[2] == speeds
        return cache[3]
    end
    check_jet_data(heights, speeds)
    coeffs = fit_jet(heights, speeds)
    s.jet_cache = (Vector{Float64}(heights), Vector{Float64}(speeds), coeffs)
    coeffs
end

@inline function calc_wind_factor6(s::AM, height); jet_model(height, jet_coeffs(s)) / s.set.v_wind; end
@inline function calc_wind_factor(s::AM, height, ::Type{Val{6}}); calc_wind_factor6(s, height); end

"""
    calc_wind_factor(am::AM, height; profile_law::Int64=am.set.profile_law)

Calculates the wind factor at a given `height` using the specified wind profile law.

# Arguments
- `am::AM`: An instance of the `AM` type containing atmospheric model parameters.
- `height`: The height (in meters) at which to calculate the wind factor.
- `profile_law::Int64`: (Optional) The wind profile law to use for the calculation. 
  Defaults to `am.set.profile_law`.

# Returns
- The wind factor at the specified height as determined by the chosen profile law.
"""
@inline function calc_wind_factor(am::AM, height, profile_law::Int64=am.set.profile_law)
    if profile_law == 0
        1.0
    elseif profile_law == 1
        calc_wind_factor1(am, height)
    elseif profile_law == 2
        calc_wind_factor2(am, height)
    elseif profile_law == 3
        calc_wind_factor3(am, height)
    elseif profile_law == 4
        calc_wind_factor4(am, height)
    elseif profile_law == 5
        calc_wind_factor5(am, height)
    elseif profile_law == 6
        calc_wind_factor6(am, height)
    else
        throw(DomainError(profile_law, "invalid profile_law"))
    end
end

include("windfield.jl")

end