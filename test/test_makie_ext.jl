using MakieControlPlots
using MakieControlPlots.Makie: Contour, Colorbar

"True if `axis[1:n]` spans at most `height` and `axis[1:n+1]` would span more."
cropped_to(axis, n, height) = axis[n] - first(axis) <= height &&
                              (n == length(axis) || axis[n + 1] - first(axis) > height)

@testset "plot(am) shows ±2σ isosurfaces of the scaled along-wind turbulence" begin
    fig = plot(am)
    gusts = only(p for p in fig.content[1].scene.plots if p isa Contour)
    turbulence = gusts[4][]
    x, y, z = AtmosphericModels.grid_axes(am)
    nx, ny = size(turbulence)
    height = last(z) - first(z)
    @test cropped_to(x, nx, height)
    @test cropped_to(y, ny, height)
    @test collect(gusts[1][]) ≈ [first(x), x[nx]]
    @test collect(gusts[2][]) ≈ [first(y), y[ny]]
    @test collect(gusts[3][]) ≈ [first(z), last(z)]
    scale = am.set.use_turbulence * rel_turbo(am, am.wf.v_wind_gnd)
    @test AtmosphericModels.turbulence_scale(am, am.wf) == scale
    u = am.wf.u[1:nx, 1:ny, :] .* scale
    @test turbulence ≈ u
    σ = sqrt(sum(abs2, u) / length(u))
    @test gusts.levels[] ≈ [-2σ, 2σ]
    @test 0 < gusts.alpha[] < 1
    limit = maximum(abs, u)
    @test collect(gusts.colorrange[]) ≈ [-limit, limit]
    colorbar = only(c for c in fig.content if c isa Colorbar)
    @test collect(colorbar.limits[]) ≈ [-limit, limit]
end

@testset "plot(am) without a wind field throws" begin
    @test_throws ArgumentError plot(AtmosphericModel(set; nowindfield=true))
end
