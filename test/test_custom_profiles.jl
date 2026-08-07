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

@testset "custom_jet" begin
    heights = [10.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0]
    c_bg, a_bg, U_J, z_c, sigma = 2.0, 0.2, 5.0, 200.0, 30.0
    speeds = [c_bg * z^a_bg + U_J * exp(-(z - z_c)^2 / (2 * sigma^2)) for z in heights]

    @test custom_jet(heights, speeds, 200.0) ≈ speeds[5] rtol=1e-6
    @test custom_jet(heights, speeds, 400.0) ≈ speeds[8] rtol=1e-6

    v200 = custom_jet(heights, speeds, 200.0)
    v250 = custom_jet(heights, speeds, 250.0)
    v300 = custom_jet(heights, speeds, 300.0)
    @test v200 > v250 > v300

    @test_throws ArgumentError custom_jet(heights, speeds[1:2], 200.0)
    @test_throws ArgumentError custom_jet(heights[1:4], speeds[1:4], 200.0)
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

    heights_jet = [10.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0]
    c_bg, a_bg, U_J, z_c, sigma = 2.0, 0.2, 5.0, 200.0, 30.0
    speeds_jet = [c_bg * z^a_bg + U_J * exp(-(z - z_c)^2 / (2 * sigma^2)) for z in heights_jet]
    am.set.heights = heights_jet
    am.set.speeds = speeds_jet
    am.set.v_wind = speeds_jet[1]
    @test calc_wind_factor(am, 200.0, Int(CUSTOM_JET)) ≈
        custom_jet(heights_jet, speeds_jet, 200.0) / am.set.v_wind
    @test calc_wind_factor(am, 200.0, Val{Int(CUSTOM_JET)}) ≈
        calc_wind_factor(am, 200.0, Int(CUSTOM_JET))

    am.set.v_wind = v_wind_orig
    am.set.heights = [6.0]
    am.set.speeds = [am.set.v_wind]
end
nothing
