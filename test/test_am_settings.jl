environment_file = joinpath(@__DIR__, "..", "data", "settings.yaml")

function same_fields(amset::AMSettings, set::Settings)
    all(getproperty(amset, field) == getproperty(set, field)
        for field in fieldnames(AMSettings))
end

@testset "AMSettings" begin
    @testset "defaults match KiteUtils.Settings" begin
        @test same_fields(AMSettings(), Settings())
    end

    reference = Settings("system.yaml"; relax=true)

    @testset "reads the environment block that load_settings reads" begin
        amset = AMSettings(environment_file)
        @test same_fields(amset, reference)
        @test amset.grid isa Vector{Int64}
    end

    @testset "full KiteUtils settings.yaml loads as KiteUtils loads it" begin
        kiteutils_data = joinpath(pkgdir(KiteUtils), "data")
        full_set = Settings(joinpath(kiteutils_data, "system.yaml"))
        @test same_fields(AMSettings(joinpath(kiteutils_data, "settings.yaml")), full_set)
    end

    @testset "environment-only yaml: wind_vec and profile_law as in load_settings" begin
        mktempdir() do dir
            file = joinpath(dir, "atmosphere_settings.yaml")
            write(file, """
                environment:
                    v_wind: 9.0
                    upwind_dir: -90.0
                    wind_vec: [3.0, -4.0, 0.0]
                    use_wind_vec: true
                    profile_law: 1
                """)
            expected = Settings()
            expected.use_wind_vec = true
            expected.wind_vec = [3.0, -4.0, 0.0]
            amset = AMSettings(file)
            @test amset.v_wind ≈ expected.v_wind
            @test amset.upwind_dir ≈ expected.upwind_dir
            @test amset.profile_law == 1

            write(file, "environment:\n    profile_law: 7\n")
            @test_throws ArgumentError AMSettings(file)
        end
    end

    @testset "AtmosphericModel from AMSettings matches the one from Settings" begin
        am_ref = AtmosphericModel(reference)
        am_env = AtmosphericModel(AMSettings(environment_file))
        @test am_env.set isa AMSettings
        @test am_env.rho_zero_temp == am_ref.rho_zero_temp
        @test calc_rho(am_env, 100.0) == calc_rho(am_ref, 100.0)
        for law in (EXP, LOG, EXPLOG)
            @test calc_wind_factor(am_env, 100.0, Int(law)) ==
                  calc_wind_factor(am_ref, 100.0, Int(law))
        end
        position_time = (20.0, 10.0, 150.0, 3.0)
        @test get_wind(am_env, position_time...) == get_wind(am_ref, position_time...)
        @test AtmosphericModels.calc_basename(am_env.set) ==
              AtmosphericModels.calc_basename(reference)
    end

    @testset "AtmosphericModel(set::Settings) keeps sharing the caller's Settings" begin
        am_shared = AtmosphericModel(reference; nowindfield=true)
        @test am_shared.set === reference
    end
end
