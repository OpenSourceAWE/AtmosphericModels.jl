# set_data_path("data") 

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
    datapath = get_data_path()
    tmpdir = joinpath(mktempdir(cleanup=true), "data")
    mkpath(tmpdir)
    olddir = pwd()
    cd(dirname(tmpdir))
    set_data_path(tmpdir)

    v_wind_gnd = 5.324
    fullname = AtmosphericModels.calcFullName(v_wind_gnd)
    @test basename(fullname) == "windfield_4050_500_1.0_5.3"
    x, y, z = create_xyz()
    u, v, w = create_uvw()
    param = [1, 2]
    AtmosphericModels.save(x, y, z, u, v, w, param; v_wind_gnd)
    println("Saved windfield data to: ", fullname*".npz")
    @test isfile(fullname * ".npz")

    # Load the data back
    x2, y2, z2, u2, v2, w2, param2 = AtmosphericModels.load(;v_wind_gnd)
    @test x == x2
    @test y == y2
    @test z == z2
    @test u ≈ u2
    @test v ≈ v2
    @test w ≈ w2
    @test param == param2

    am = AtmosphericModel(set)
    windfield = AtmosphericModels.load_windfield(am, v_wind_gnd+0.2)
    @test typeof(windfield) == Tuple{Vector{Int64}, Vector{Int64}, Vector{Int64}, Array{Float64, 3}, Array{Float64, 3}, Array{Float64, 3}, Vector{Int64}}

    grid = AtmosphericModels.create_grid()
    @test typeof(grid) == Tuple{Array{Float64, 3}, Array{Float64, 3}, Array{Float64, 3}}

    x = range(0, 50, length=25)
    y = range(0, 800, length=400)
    z = range(0, 200, length=100)

    u, v, w = AtmosphericModels.create_windfield(x, y, z; sigma1=1.2)
    am = AtmosphericModel(set)
    am.set.v_wind = v_wind_gnd
    AtmosphericModels.addWindSpeed(am, z, u)

    set_data_path(olddir)
    cd(olddir)
end