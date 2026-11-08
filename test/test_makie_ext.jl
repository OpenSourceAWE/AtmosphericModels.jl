using MakieControlPlots
using MakieControlPlots: Makie
using MakieControlPlots.Makie: Axis3, Contour, SliderGrid, set_close_to!

"True if `axis[1:n]` spans at most `height` and `axis[1:n+1]` would span more."
cropped_to(axis, n, height) = axis[n] - first(axis) <= height &&
                              (n == length(axis) || axis[n + 1] - first(axis) > height)

"The gust and lull isosurfaces of the figure `fig`."
isosurfaces(fig) =
    [p for p in only(c for c in fig.content if c isa Axis3).scene.plots if p isa Contour]

"Scaled along-wind turbulence plotted by `plot(am)` and its standard deviation."
function plotted_turbulence(am, nx, ny)
    u = am.wf.u[1:nx, 1:ny, :] .* AtmosphericModels.turbulence_scale(am, am.wf)
    return u, sqrt(sum(abs2, u) / length(u))
end

@testset "plot(am) shows red gusts and blue lulls beyond ±2σ of the scaled u'" begin
    gust, lull = isosurfaces(plot(am))
    turbulence = gust[4][]
    x, y, z = AtmosphericModels.grid_axes(am)
    nx, ny = size(turbulence)[1:2]
    height = last(z) - first(z)
    @test cropped_to(x, nx, height)
    @test cropped_to(y, ny, height)
    @test collect(gust[1][]) ≈ [first(x), x[nx]]
    @test collect(gust[2][]) ≈ [first(y), y[ny]]
    @test collect(gust[3][]) ≈ [first(z), last(z)]
    scale = am.set.use_turbulence * rel_turbo(am, am.wf.v_wind_gnd)
    @test AtmosphericModels.turbulence_scale(am, am.wf) == scale
    u, σ = plotted_turbulence(am, nx, ny)
    @test turbulence ≈ u
    @test lull[4][] ≈ u
    @test gust.levels[] ≈ [2σ]
    @test lull.levels[] ≈ [-2σ]
    @test all(==(Makie.to_color(:red)), Makie.to_colormap(gust.colormap[]))
    @test all(==(Makie.to_color(:blue)), Makie.to_colormap(lull.colormap[]))
    @test 0 < gust.alpha[] < 1
end

@testset "plot(am; threshold) puts the isosurfaces at ±threshold" begin
    gust, lull = isosurfaces(plot(am; threshold=0.7))
    @test gust.levels[] == [0.7]
    @test lull.levels[] == [-0.7]
end

@testset "plot_interactive(am) moves the isosurfaces with the threshold slider" begin
    fig = plot_interactive(am)
    gust, lull = isosurfaces(fig)
    u, σ = plotted_turbulence(am, size(gust[4][])[1:2]...)
    slider = only(only(c for c in fig.content if c isa SliderGrid).sliders)
    @test slider.value[] ≈ 2σ atol=step(slider.range[])
    @test last(slider.range[]) ≈ maximum(abs, u)
    set_close_to!(slider, 0.5)
    @test gust.levels[] == [slider.value[]]
    @test lull.levels[] == [-slider.value[]]
    @test slider.value[] ≈ 0.5 atol=step(slider.range[])
end

@testset "plot(am) without a wind field throws" begin
    @test_throws ArgumentError plot(AtmosphericModel(set; nowindfield=true))
end
