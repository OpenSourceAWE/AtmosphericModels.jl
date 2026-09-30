module AtmosphericModelsMakieExt

using MakieControlPlots.Makie
import MakieControlPlots
using AtmosphericModels
using AtmosphericModels: grid_axes, turbulence_scale

"Isosurface level of `plot(am)`, in standard deviations of u'."
const ISO_SIGMAS = 2
"Opacity of the isosurfaces of `plot(am)`."
const ISO_ALPHA = 0.5

"""
    plot(am::AtmosphericModel)

Plot the turbulent wind field of `am` in 3D as see-through isosurfaces of the along-wind
turbulence u' [m/s], scaled as `get_wind` applies it, at ±2 standard deviations: the
gusts and lulls inside the grid. The long horizontal axis is cut to the height of the
grid. Returns the `Figure`; needs a field in `am.wf`.
"""
function MakieControlPlots.plot(am::AtmosphericModel)
    wf = am.wf
    wf === nothing && throw(ArgumentError(
        "no wind field: AtmosphericModel(set) only loads one when set.use_turbulence > 0"))
    x, y, z = grid_axes(am)
    height = last(z) - first(z)
    nx = findlast(<=(first(x) + height), x)
    ny = findlast(<=(first(y) + height), y)
    x, y = x[1:nx], y[1:ny]
    turbulence = wf.u[1:nx, 1:ny, :] .* turbulence_scale(am, wf)
    limit = maximum(abs, turbulence)
    level = ISO_SIGMAS * sqrt(sum(abs2, turbulence) / length(turbulence))

    fig = Figure(; size=(900, 650))
    ax = Axis3(fig[1, 1]; aspect=:data, azimuth=-0.3π, elevation=0.2π,
               xlabel="x [m]", ylabel="y [m]", zlabel="z [m]",
               title="Turbulent wind field, v_wind_gnd = $(wf.v_wind_gnd) m/s")
    contour!(ax, extrema(x), extrema(y), extrema(z), turbulence; levels=[-level, level],
             colormap=:balance, colorrange=(-limit, limit), alpha=ISO_ALPHA)
    Colorbar(fig[1, 2]; colormap=:balance, limits=(-limit, limit),
             label="along-wind turbulence u' [m/s]")
    return fig
end

end
