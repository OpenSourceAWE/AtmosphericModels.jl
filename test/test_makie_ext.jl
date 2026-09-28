using MakieControlPlots
using MakieControlPlots.Makie: Surface, Colorbar

@testset "plot(am) shows the scaled along-wind turbulence on three faces" begin
    fig = plot(am)
    faces = [p for p in fig.content[1].scene.plots if p isa Surface]
    @test length(faces) == 3
    x, y, z = AtmosphericModels.grid_axes(am)
    rel_turb = am.set.use_turbulence * rel_turbo(am, am.wf.v_wind_gnd)
    nx = findlast(<=(last(z) - first(z)), x)
    ny = length(y)
    @test faces[1].color[] ≈ am.wf.u[nx, 1:ny, :] .* rel_turb
    @test faces[2].color[] ≈ am.wf.u[1:nx, 1, :] .* rel_turb
    @test faces[3].color[] ≈ am.wf.u[1:nx, 1:ny, end] .* rel_turb
    limit = maximum(abs, am.wf.u[1:nx, 1:ny, :]) * rel_turb
    @test all(face -> collect(face.colorrange[]) ≈ [-limit, limit], faces)
    @test count(c -> c isa Colorbar, fig.content) == 1
end

@testset "plot(am) without a wind field throws" begin
    @test_throws ArgumentError plot(AtmosphericModel(set; nowindfield=true))
end
