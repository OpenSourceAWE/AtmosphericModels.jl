using AtmosphericModels, KiteUtils, BenchmarkTools
using Test

cd("..")
KiteUtils.set_data_path() 
set = load_settings("system.yaml")
am = AtmosphericModel(set)

include("test_windfield.jl")

@testset "calc_wind_factor" begin
    @test calc_wind_factor(am, 6.0, Val{Int(EXP)}) ≈ 1.0
    @test calc_wind_factor(am, 6.0, Val{Int(LOG)}) ≈ 1.0
    @test calc_wind_factor(am, 6.0, Val{Int(EXPLOG)}) ≈ 1.0
    @test abs(calc_wind_factor(am, 6.0, Val{Int(FAST_EXP)}) - 1.0) < 1.5e-5
    @test abs(calc_wind_factor(am, 6.0, Val{Int(FAST_LOG)}) - 1.0) < 1.5e-5
    @test abs(calc_wind_factor(am, 6.0, Val{Int(FAST_EXPLOG)}) - 1.0) < 1.5e-5    
    heights = 6:10:1000
    for height in heights
        ref = calc_wind_factor(am, height, Val{Int(1)})
        approx = calc_wind_factor(am, height, Val{Int(4)})
        @test abs(approx/ref - 1.0) < 1.5e-5
    end
    for height in heights
        ref = calc_wind_factor(am, height, Val{Int(2)})
        approx = calc_wind_factor(am, height, Val{Int(5)})
        @test abs(approx/ref - 1.0) < 1.5e-5
    end
    for height in heights
        ref = calc_wind_factor(am, height, Val{Int(3)})
        approx = calc_wind_factor(am, height, Val{Int(6)})
        @test abs(approx/ref - 1.0) < 1.5e-5
    end
end

@testset "calc_rho        " begin
    @test calc_rho(am, 0.0) ≈ am.set.rho_0
    am.set.temp_ref = 15 - AtmosphericModels.ABS_ZERO + 15.0
    clear(am)
    @test calc_rho(am, 0.0) ≈ 0.5 * am.set.rho_0
    am.set.temp_ref = 15
    clear(am)
end


