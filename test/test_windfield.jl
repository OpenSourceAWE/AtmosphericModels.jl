# set_data_path("data") 

@testset "windfield       " begin
    @test AtmosphericModels.pfq(0.1) ≈ 1.079576584249971
    @test AtmosphericModels.calc_sigma1(am, 10.0) ≈ 3.1692995457170285
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

@testset "3d_windfield    " begin
    datapath = get_data_path()
    tmpdir = joinpath(mktempdir(cleanup=true), "data")
    mkpath(tmpdir)
    olddir = pwd()
    cd(dirname(tmpdir))
    set_data_path(tmpdir)

    v_wind_gnd = 5.324
    fullname = AtmosphericModels.calc_full_name(v_wind_gnd)
    @test basename(fullname) == "windfield_4050_500_1.0_5.3"
    x, y, z = create_xyz()
    u, v, w = create_uvw()
    param = [1, 2]
    am = AtmosphericModel(set)
    AtmosphericModels.save(am, x, y, z, u, v, w, param; v_wind_gnd)
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

    y, x, z = AtmosphericModels.create_grid(10, 20, 10, 5)

    u, v, w = AtmosphericModels.create_windfield(x, y, z; sigma1=1.2)
    am = AtmosphericModel(set)
    am.set.v_wind = v_wind_gnd
    AtmosphericModels.addWindSpeed(am, z, u)

    set_data_path(olddir)
    cd(olddir)
end

y, x, z = AtmosphericModels.create_grid(10, 20, 10, 5)
u, v, w = AtmosphericModels.create_windfield(x, y, z, sigma1=1.0)

@testset "create_windfield" begin
    nx, ny, nz = size(x)
    @test nx == 11
    @test ny == 6
    @test nz == 6
    # Domain lengths
    Lx = x[end,1,1] - x[1,1,1]
    Ly = y[1,end,1] - y[1,1,1]
    Lz = z[1,1,end] - z[1,1,1]
    @test Lx == 20.0
    @test Ly == 10.0
    @test Lz == 10.0
    @test AtmosphericModels.pfq(0.5) ≈ 1.7936563627777333
    @test sum(x) == 3960.0
    @test sum(y) == 0.0
    @test sum(z) == 3960.0
    @test all(x[1,:,:] .== 0)
    @test all(x[2,:,:] .== 2)
    @test all(x[11,:,:] .== 20)
    @test size(u) == (11,6,6)
    @test size(v) == (11,6,6)
    @test size(w) == (11,6,6)
    @test std(u) ≈ 1.0
    @test std(v) ≈ 0.7
    @test std(w) ≈ 0.5
end