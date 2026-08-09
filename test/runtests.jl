using AtmosphericModels, KiteUtils, BenchmarkTools
using Test

# Absolute, so that the suite can be included from anywhere; nothing here depends on the
# working directory.
KiteUtils.set_data_path(normpath(joinpath(@__DIR__, "..", "data")))
set = load_settings("system.yaml"; relax=true)
am = AtmosphericModel(set)

include("test_windfield.jl")
include("test_custom_profiles.jl")

@testset "calc_wind_factor" begin
    @test calc_wind_factor(am, 6.0, Int(CONSTANT)) == 1.0
    @test calc_wind_factor(am, 100.0, Int(CONSTANT)) == 1.0
    @test calc_wind_factor(am, 100.0, Val{Int(CONSTANT)}) == 1.0
    @test calc_wind_factor(am, 6.0, Val{Int(EXP)}) ≈ 1.0
    @test calc_wind_factor(am, 6.0, Val{Int(LOG)}) ≈ 1.0
    @test calc_wind_factor(am, 6.0, Val{Int(EXPLOG)}) ≈ 1.0
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


