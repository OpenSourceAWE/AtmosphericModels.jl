@testset "windfield" begin
    @test AtmosphericModels.pfq(0.1) ≈ 1.079576584249971
    @test AtmosphericModels.calc_sigma1(am, 10.0) ≈ 2.181983002542761
    @test AtmosphericModels.nextpow2(10) == 16
end

@testset "3d_windfield" begin
    fullname = AtmosphericModels.calcFullName(6.0)
    @test basename(fullname) == "windfield_4050_500_1.0_6.0"
end