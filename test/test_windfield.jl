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
    fullname = AtmosphericModels.calc_full_name(v_wind_gnd; basename="windfield_4050_100_500_70")
    @test basename(fullname) == "windfield_4050_100_500_70_1.0_5.3"
    x, y, z = create_xyz()
    u, v, w = create_uvw()
    param = [1, 2]
    am = AtmosphericModel(set; nowindfield=true)
    AtmosphericModels.save(am, x, y, z, u, v, w, param; v_wind_gnd=am.set.v_wind)
    @test isfile(fullname * ".npz")

    # Load the data back
    x2, y2, z2, u2, v2, w2, param2 = AtmosphericModels.load(am;v_wind_gnd)
    @test x == x2
    @test y == y2
    @test z == z2
    @test u ≈ u2
    @test v ≈ v2
    @test w ≈ w2
    @test param == param2

    am = AtmosphericModel(set; nowindfield=true)
    windfield = AtmosphericModels.load_windfield(am, v_wind_gnd+0.2)
    @test typeof(windfield) == Tuple{Vector{Int64}, Vector{Int64}, Vector{Int64}, Array{Float64, 3}, Array{Float64, 3},
                                     Array{Float64, 3}, Vector{Int64}}

    am.set.grid=[20, 10, 10, 5]
    grid = AtmosphericModels.create_grid(am)
    @test typeof(grid) == Tuple{Array{Float64, 3}, Array{Float64, 3}, Array{Float64, 3}}

    y, x, z = AtmosphericModels.create_grid(am)
    u, v, w = AtmosphericModels.create_windfield(x, y, z; sigma1=1.2)

    set_data_path(olddir)
    cd(olddir)
end

am.set.grid=[20, 10, 10, 5]
y, x, z = AtmosphericModels.create_grid(am)
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

@testset "calc_turbulent_wind" begin
    # Reference values pinned against the pre-move KiteModels.calc_turbulent_wind
    # implementation (same rotation, Taylor advection and nearest-grid lookup),
    # captured before that implementation was deleted.
    v_wind, v_wind_tether = AtmosphericModels.calc_turbulent_wind(am, [10.0, 20.0, 100.0], 12.5; upwind_dir=0.3)
    @test v_wind ≈ [-2.144770095983673, -11.124313263469531, -0.6070299938532554]
    @test v_wind_tether ≈ [-1.4502831306377473, -8.76678510107124, -0.3167369264079984]

    v_wind2, v_wind_tether2 = AtmosphericModels.calc_turbulent_wind(am, [0.0, 0.0, 50.0], 0.0; upwind_dir=-pi/2)
    @test v_wind2 ≈ [9.796027271006283, -0.34202755189268585, -0.3957750583820598]
    @test v_wind_tether2 ≈ [8.678228770857402, 0.40297265466661003, -0.07675347144301002]

    # use_turbulence == 0: deterministic mean wind, no wind field needed
    set0 = deepcopy(am.set)
    set0.use_turbulence = 0.0
    am0 = AtmosphericModel(set0)
    v_wind0, v_wind_tether0 = AtmosphericModels.calc_turbulent_wind(am0, [0.0, 0.0, 100.0], 5.0; upwind_dir=0.0)
    v_height = set0.v_wind * calc_wind_factor(am0, 100.0)
    v_height_half = set0.v_wind * calc_wind_factor(am0, 50.0)
    @test v_wind0[1] ≈ 0.0 atol=1e-10
    @test v_wind0[2] ≈ -v_height
    @test v_wind0[3] == 0.0
    @test v_wind_tether0[1] ≈ 0.0 atol=1e-10
    @test v_wind_tether0[2] ≈ -v_height_half
    @test v_wind_tether0[3] == 0.0
end