@testset "custom_log/custom_exp" begin
    heights = [10.0, 50.0, 100.0, 200.0, 400.0]

    a_log, b_log = 1.5, 2.0
    speeds_log = a_log .* log.(heights) .+ b_log
    @test custom_log(heights, speeds_log, 150.0) ≈ a_log * log(150.0) + b_log
    @test custom_log(heights, speeds_log, heights[3]) ≈ speeds_log[3]

    c_exp, a_exp = 2.0, 0.3
    speeds_exp = c_exp .* heights .^ a_exp
    @test custom_exp(heights, speeds_exp, 150.0) ≈ c_exp * 150.0^a_exp
    @test custom_exp(heights, speeds_exp, heights[3]) ≈ speeds_exp[3]

    @test_throws ArgumentError custom_log(heights, speeds_log[1:2], 150.0)
    @test_throws ArgumentError custom_log(heights[1:1], speeds_log[1:1], 150.0)
    @test_throws ArgumentError custom_exp(heights, speeds_exp[1:2], 150.0)
    @test_throws ArgumentError custom_exp(heights[1:1], speeds_exp[1:1], 150.0)
end

@testset "calc_wind_factor CUSTOM_*" begin
    v_wind_orig = am.set.v_wind
    heights = [10.0, 50.0, 100.0, 200.0, 400.0]
    a_log, b_log = 1.5, 2.0
    speeds_log = a_log .* log.(heights) .+ b_log
    am.set.heights = heights
    am.set.speeds = speeds_log
    am.set.v_wind = speeds_log[1]
    @test calc_wind_factor(am, 150.0, Int(CUSTOM_LOG)) ≈
        custom_log(heights, speeds_log, 150.0) / am.set.v_wind
    @test calc_wind_factor(am, 150.0, Val{Int(CUSTOM_LOG)}) ≈
        calc_wind_factor(am, 150.0, Int(CUSTOM_LOG))

    c_exp, a_exp = 2.0, 0.3
    speeds_exp = c_exp .* heights .^ a_exp
    am.set.speeds = speeds_exp
    am.set.v_wind = speeds_exp[1]
    @test calc_wind_factor(am, 150.0, Int(CUSTOM_EXP)) ≈
        custom_exp(heights, speeds_exp, 150.0) / am.set.v_wind
    @test calc_wind_factor(am, 150.0, Val{Int(CUSTOM_EXP)}) ≈
        calc_wind_factor(am, 150.0, Int(CUSTOM_EXP))

    @test_throws ErrorException calc_wind_factor(am, 150.0, Int(CUSTOM_JET))

    am.set.v_wind = v_wind_orig
    am.set.heights = [6.0]
    am.set.speeds = [am.set.v_wind]
end
