module AtmosphericModels

using KiteUtils
using HypergeometricFunctions:_₂F₁
using NPZ, Printf
using FFTW, LinearAlgebra, Random, Statistics

export AtmosphericModel, ProfileLaw, EXP, LOG, EXPLOG
export clear, calc_rho, calc_wind_factor, rel_turbo

export new_windfield, new_windfields, get_wind

const ABS_ZERO = -273.15
const SRL = StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}

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
"""
Base.@kwdef mutable struct AtmosphericModel
    set::Settings
    rho_zero_temp::Float64 = (15.0 - ABS_ZERO) / (set.temp_ref - ABS_ZERO) * set.rho_0
    wf::Union{WindField, Nothing} = nothing
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
    @enum ProfileLaw EXP=1 LOG=2 EXPLOG=3

Enumeration to describe the wind profile low that is used.
""" ProfileLaw

@enum ProfileLaw begin 
    EXP=1 
    LOG=2 
    EXPLOG=3 
end

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


# Calculate the wind speed at a given height and reference height.
@inline function calc_wind_factor1(s::AM, height);  exp(s.set.alpha * log(height/s.set.h_ref)); end
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
    if profile_law == 1
        calc_wind_factor1(am, height)
    elseif profile_law == 2
        calc_wind_factor2(am, height)
    elseif profile_law == 3
        calc_wind_factor3(am, height)   
    else
        throw(DomainError(profile_law, "invalid profile_law"))
    end
end

include("windfield.jl")

end