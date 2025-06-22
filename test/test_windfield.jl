@testset "windfield" begin
    @test AtmosphericModels.pfq(0.1) ≈ 1.079576584249971
    @test AtmosphericModels.calc_sigma1(am, 10.0) ≈ 2.181983002542761
    @test AtmosphericModels.nextpow2(10) == 16
end

function create_xyz()
    x = 0:10:100
    y = 0:10:200
    z = 0:10:300
    return x, y, z
end

function create_uvw()
    u = rand(11, 21, 31)
    v = rand(11, 21, 31)
    w = rand(11, 21, 31)
    return u, v, w
end

@testset "3d_windfield" begin
    fullname = AtmosphericModels.calcFullName(6.0)
    @test basename(fullname) == "windfield_4050_500_1.0_6.0"
    x, y, z = create_xyz()
    u, v, w = create_uvw()
    param = [1, 2]
    AtmosphericModels.save(x, y, z, u, v, w, param)
end