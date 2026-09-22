"""
    mutable struct AMSettings

The `environment:` settings of `KiteUtils.Settings` that an [`AtmosphericModel`](@ref)
reads, with the same names, units and defaults.
"""
Base.@kwdef mutable struct AMSettings
    "wind speed at reference height [m/s]"
    v_wind::Float64 = 0.0
    "initial upwind direction [deg]"
    upwind_dir::Float64 = 0.0
    "temperature at reference height [°C]"
    temp_ref::Float64 = 0.0
    "height of the ground station above sea level [m]"
    height_gnd::Float64 = 0.0
    "reference height for the wind speed [m]"
    h_ref::Float64 = 0.0
    "air density at zero height and 15 °C [kg/m³]"
    rho_0::Float64 = 0.0
    "exponent of the wind profile law"
    alpha::Float64 = 0.0
    "surface roughness [m]"
    z0::Float64 = 0.0
    "0=CONST, 1=EXP, 2=LOG, 3=EXPLOG, 4=CUSTOM_LOG, 5=CUSTOM_EXP, 6=CUSTOM_JET"
    profile_law::Int64 = 0
    "heights at which the wind speed is given, for the CUSTOM_* profile laws [m]"
    heights::Vector{Float64} = [6.0]
    "wind speeds at the given heights, for the CUSTOM_* profile laws [m/s]"
    speeds::Vector{Float64} = [v_wind]
    "turbulence intensity relative to Cabauw, NL"
    use_turbulence::Float64 = 0.0
    "wind speeds at ref height for calculating the turbulent wind field [m/s]"
    v_wind_gnds::Vector{Float64} = Float64[]
    "average height during reel out [m]"
    avg_height::Float64 = 0.0
    "relative turbulence at the v_wind_gnds"
    rel_turbs::Vector{Float64} = Float64[]
    "expected value of the turbulence intensity at 15 m/s"
    i_ref::Float64 = 0.0
    "five times the average wind speed at hub height over the full year [m/s]"
    v_ref::Float64 = 0.0
    "grid size nx, ny, nz and minimal height z_min [m]"
    grid::Vector{Int64} = Int64[]
    "grid resolution in z direction [m]"
    height_step::Float64 = 0.0
    "grid resolution in x and y direction [m]"
    grid_step::Float64 = 0.0
end

StructTypes.StructType(::Type{AMSettings}) = StructTypes.Mutable()

"""
    AMSettings(file::AbstractString)

Load the `environment:` section of the yaml `file`, in the schema of KiteUtils'
`settings.yaml`; other sections and keys are ignored, missing keys keep their defaults.
"""
function AMSettings(file::AbstractString)
    set = AMSettings()
    KiteUtils.update_settings(YAML.load_file(file), ["environment"], set)
    set
end

const AtmosphereSettings = Union{Settings, AMSettings}
