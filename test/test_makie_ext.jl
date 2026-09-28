using MakieControlPlots
using MakieControlPlots.Makie: Surface, Colorbar

"True if `axis[1:n]` spans at most `height` and `axis[1:n+1]` would span more."
cropped_to(axis, n, height) =
    axis[n] - first(axis) <= height && (n == length(axis) || axis[n + 1] - first(axis) > height)

@testset "plot(am) shows the scaled along-wind turbulence on three faces" begin
    fig = plot(am)
    faces = [p for p in fig.content[1].scene.plots if p isa Surface]
    @test length(faces) == 3
    x, y, z = AtmosphericModels.grid_axes(am)
    nx, ny = size(faces[3].color[])
    height = last(z) - first(z)
    @test cropped_to(x, nx, height)
    @test cropped_to(y, ny, height)
    scale = am.set.use_turbulence * rel_turbo(am, am.wf.v_wind_gnd)
    @test AtmosphericModels.turbulence_scale(am, am.wf) == scale
    @test faces[1].color[] ≈ am.wf.u[nx, 1:ny, :] .* scale
    @test faces[2].color[] ≈ am.wf.u[1:nx, 1, :] .* scale
    @test faces[3].color[] ≈ am.wf.u[1:nx, 1:ny, end] .* scale
    limit = maximum(abs, am.wf.u[1:nx, 1:ny, :]) * scale
    @test all(face -> collect(face.colorrange[]) ≈ [-limit, limit], faces)
    @test count(c -> c isa Colorbar, fig.content) == 1
end

@testset "plot(am) without a wind field throws" begin
    @test_throws ArgumentError plot(AtmosphericModel(set; nowindfield=true))
end
