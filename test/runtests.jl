using AtmosphericModels, KiteUtils, BenchmarkTools
using Test

# Absolute, so that the suite can be included from anywhere; nothing here depends on the
# working directory.
KiteUtils.set_data_path(normpath(joinpath(@__DIR__, "..", "data")))
set = load_settings("system.yaml"; relax=true)
am = AtmosphericModel(set)

include("test_windfield.jl")
include("test_custom_profiles.jl")
include("test_am_settings.jl")

@testset "calc_wind_factor" begin
    @test calc_wind_factor(am, 6.0, Int(CONSTANT)) == 1.0
    @test calc_wind_factor(am, 100.0, Int(CONSTANT)) == 1.0
    @test calc_wind_factor(am, 100.0, Val{Int(CONSTANT)}) == 1.0
    @test calc_wind_factor(am, 6.0, Val{Int(EXP)}) ≈ 1.0
    @test calc_wind_factor(am, 6.0, Val{Int(LOG)}) ≈ 1.0
    @test calc_wind_factor(am, 6.0, Val{Int(EXPLOG)}) ≈ 1.0

    # Int64-dispatched runtime overload: each branch must agree with its
    # Val-dispatched counterpart at the reference height and at another height.
    for law in (CONSTANT, EXP, LOG, EXPLOG)
        @test calc_wind_factor(am, 6.0, Int(law)) ≈ calc_wind_factor(am, 6.0, Val{Int(law)})
        @test calc_wind_factor(am, 100.0, Int(law)) ≈ calc_wind_factor(am, 100.0, Val{Int(law)})
    end

    # default profile_law argument falls back to am.set.profile_law
    @test calc_wind_factor(am, 100.0) == calc_wind_factor(am, 100.0, am.set.profile_law)

    # an out-of-range profile_law must raise a DomainError
    @test_throws DomainError calc_wind_factor(am, 6.0, -1)
    @test_throws DomainError calc_wind_factor(am, 6.0, 7)
end

@testset "calc_rho        " begin
    @test calc_rho(am, 0.0) ≈ am.set.rho_0
    am.set.temp_ref = 15 - AtmosphericModels.ABS_ZERO + 15.0
    clear(am)
    @test calc_rho(am, 0.0) ≈ 0.5 * am.set.rho_0
    am.set.temp_ref = 15
    clear(am)
end
nothing


