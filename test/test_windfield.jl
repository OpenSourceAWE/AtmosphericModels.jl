# set_data_path("data") 

@testset "windfield       " begin
    @test AtmosphericModels.pfq(0.1) ≈ 1.079576584249971
    @test AtmosphericModels.calc_sigma1(am, 10.0) ≈ 3.1692995457170285
    @test AtmosphericModels.nextpow2(10) == 16

    # The axes are rebuilt from set.grid = [4050, 100, 500, 70], not read from the file.
    @test am.wf.x == first(AtmosphericModels.grid_axes(am))
    @test am.wf.x_max == 4050.0
    @test am.wf.z_min == 70.0
    @test am.wf.y_range == 100.0
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
    # Without this the dummy field below would overwrite the real one in the scratchspace.
    set_windfield_path!(tmpdir)

    v_wind_gnd = 5.324
    fullname = AtmosphericModels.calc_full_name(v_wind_gnd; basename=AtmosphericModels.calc_basename(set))
    @test startswith(basename(fullname), "windfield_4050_100_500_70_")
    @test endswith(basename(fullname), "_5.3")
    u, v, w = create_uvw()
    param = [1, 2]
    am = AtmosphericModel(set; nowindfield=true)
    AtmosphericModels.save(am, u, v, w, param; v_wind_gnd=am.set.v_wind)
    @test isfile(fullname * ".npz")

    # Load the data back
    u2, v2, w2, param2 = AtmosphericModels.load(am;v_wind_gnd)
    @test u ≈ u2
    @test v ≈ v2
    @test w ≈ w2
    @test param == param2

    am = AtmosphericModel(set; nowindfield=true)
    windfield = AtmosphericModels.load_windfield(am, v_wind_gnd+0.2)
    @test typeof(windfield) == Tuple{Array{Float64, 3}, Array{Float64, 3}, Array{Float64, 3},
                                     Vector{Int64}, Float64}
    @test windfield[end] == 5.324
    # The coordinate meshgrids are not in the file, so they cannot be read back.
    @test !haskey(AtmosphericModels.NPZ.npzread(fullname * ".npz"), "x")

    am.set.grid=[20, 10, 10, 5]
    grid = AtmosphericModels.create_grid(am)
    @test typeof(grid) == Tuple{Array{Float64, 3}, Array{Float64, 3}, Array{Float64, 3}}

    y, x, z = AtmosphericModels.create_grid(am)
    u, v, w = AtmosphericModels.create_windfield(x, y, z; sigma1=1.2)

    x_range, y_range, z_range = AtmosphericModels.grid_axes(am)
    @test x[:, 1, 1] == collect(x_range)
    @test y[1, :, 1] == collect(y_range)
    @test z[1, 1, :] == collect(z_range)

    set_data_path(olddir)
    set_windfield_path!("")
    cd(olddir)
end

@testset "windfield_path  " begin
    @test isdir(windfield_path())
    @test basename(windfield_path()) == "windfields"

    set2 = deepcopy(set)
    set2.grid = [99, 98, 97, 96]
    grid_only = AtmosphericModels.grid_basename(set2)
    hashed = AtmosphericModels.calc_basename(set2)
    @test grid_only == "windfield_99_98_97_96"
    @test hashed == grid_only * "_" * AtmosphericModels.param_digest(set2)

    tmpdir = mktempdir(cleanup=true)
    olddatapath = get_data_path()
    set_data_path(tmpdir)
    try
        @test isnothing(AtmosphericModels.find_windfield(set2, 5.324))
        write(joinpath(tmpdir, grid_only * "_1.0_5.3.npz"), "")
        @test AtmosphericModels.find_windfield(set2, 5.324) == joinpath(tmpdir, grid_only * "_1.0_5.3")
        write(joinpath(tmpdir, grid_only * "_5.3.npz"), "")
        @test AtmosphericModels.find_windfield(set2, 5.324) == joinpath(tmpdir, grid_only * "_5.3")
        # The name carrying the digest wins over both older ones.
        write(joinpath(tmpdir, hashed * "_5.3.npz"), "")
        @test AtmosphericModels.find_windfield(set2, 5.324) == joinpath(tmpdir, hashed * "_5.3")
    finally
        set_data_path(olddatapath)
    end
end

@testset "settings check  " begin
    check = AtmosphericModels.check_windfield_settings
    @test isnothing(check(set))

    # The 2026-08-08 case: a settings YAML that never declared `grid`.
    set2 = deepcopy(set)
    set2.grid = Int64[]
    err = try check(set2) catch e; e end
    @test err isa ArgumentError
    @test occursin("environment.grid", err.msg)
    # ... and it is raised where the field is asked for, not later inside get_wind.
    am2 = AtmosphericModel(set2; nowindfield=true)
    @test_throws ArgumentError WindField(am2, 5.324)
    @test_throws ArgumentError AtmosphericModel(set2)

    set2 = deepcopy(set); set2.grid = [4051, 100, 500, 70]
    @test_throws ArgumentError check(set2)
    set2 = deepcopy(set); set2.height_step = 0
    @test_throws ArgumentError check(set2)
    set2 = deepcopy(set); set2.v_wind_gnds = Float64[]
    @test_throws ArgumentError check(set2)
    set2 = deepcopy(set); set2.rel_turbs = [0.342, 0.465]
    @test_throws ArgumentError check(set2)
    set2 = deepcopy(set); set2.h_ref = 0
    @test_throws ArgumentError check(set2)
end

@testset "param_digest    " begin
    set2 = deepcopy(set)
    digest = AtmosphericModels.param_digest(set2)
    @test length(digest) == 8
    for field in (:grid_step, :height_step, :i_ref, :alpha, :avg_height, :h_ref)
        old = getfield(set2, field)
        setfield!(set2, field, old * 1.1)
        @test AtmosphericModels.param_digest(set2) != digest
        setfield!(set2, field, old)
    end
    @test AtmosphericModels.param_digest(set2) == digest
    # The grid is in the name itself, not in the digest.
    set2.grid = [99, 98, 97, 96]
    @test AtmosphericModels.param_digest(set2) == digest
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
    # Reference values depend on `am`'s wind field, generated with a StableRNG(1234) seed so they stay identical across Julia versions.
    v_wind, v_wind_tether = AtmosphericModels.calc_turbulent_wind(am, [10.0, 20.0, 100.0], 12.5; upwind_dir=0.3)
    @test v_wind ≈ [-2.727005993324385, -9.252554084036728, -0.5271022045469396]
    @test v_wind_tether ≈ [-1.8852836939363489, -8.086856626440875, -0.0910792704705225]

    v_wind2, v_wind_tether2 = AtmosphericModels.calc_turbulent_wind(am, [0.0, 0.0, 50.0], 0.0; upwind_dir=-pi/2)
    @test v_wind2 ≈ [9.605009115623051, 0.14230156648318207, -0.19371077258749492]
    @test v_wind_tether2 ≈ [8.959240704788714, 0.27915918858840116, -0.41309428536554527]

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

@testset "rel_turb source " begin
    # Sampled at the origin at t = 0, so changing v_wind cannot move the sampled grid point.
    turbulence(am) = collect(AtmosphericModels.get_wind(am, 0.0, 0.0, 100.0, 0.0)) .-
                     [am.set.v_wind * calc_wind_factor(am, 100.0, am.set.profile_law), 0.0, 0.0]

    @test am.wf.v_wind_gnd == 5.324
    turb = turbulence(am)
    v_wind = am.set.v_wind
    am.set.v_wind = 8.163
    # The correction belongs to the loaded field, not to whatever v_wind the settings now hold.
    @test turbulence(am) ≈ turb
    am.set.v_wind = v_wind
    @test rel_turbo(am, am.wf.v_wind_gnd) == 0.465
end

@testset "get_wind vector " begin
    # The vector method only hoists the position-independent work out of the loop, so it must
    # return exactly what the scalar method does, not just approximately.
    positions = [KiteUtils.SVec3(10.0 + 3i, 20.0 - 2i, 50.0 + i) for i in 1:20]
    t = 12.5
    ref = [KiteUtils.SVec3(get_wind(am, p[1], p[2], p[3], t; upwind_dir=0.3)) for p in positions]

    res = get_wind(am, positions, t; upwind_dir=0.3)
    @test res isa Vector{KiteUtils.SVec3}
    @test res == ref

    # default upwind_dir, and positions that are plain vectors rather than SVec3
    @test get_wind(am, [collect(p) for p in positions], t) ==
          [KiteUtils.SVec3(get_wind(am, p[1], p[2], p[3], t)) for p in positions]

    buf = Vector{KiteUtils.SVec3}(undef, length(positions))
    @test get_wind!(buf, am, positions, t; upwind_dir=0.3) === buf
    @test buf == ref
    @test_throws DimensionMismatch get_wind!(buf, am, positions[1:end-1], t)
    @test_throws AssertionError get_wind(am, [KiteUtils.SVec3(0.0, 0.0, 4.0)], t)
end
nothing