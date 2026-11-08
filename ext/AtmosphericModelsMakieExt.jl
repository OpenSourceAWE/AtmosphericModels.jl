module AtmosphericModelsMakieExt

using MakieControlPlots.Makie
import MakieControlPlots
using AtmosphericModels
using AtmosphericModels: grid_axes, turbulence_scale

"Default threshold of `plot(am)`, in standard deviations of u'."
const THRESHOLD_SIGMAS = 2
"Opacity of the isosurfaces of `plot(am)`."
const ISO_ALPHA = 0.5
"Number of positions of the threshold slider of `plot_interactive`."
const SLIDER_STEPS = 100

"""
    plot(am::AtmosphericModel; threshold)

Plot the turbulent wind field of `am` in 3D as see-through isosurfaces of the along-wind
turbulence u' [m/s], scaled as `get_wind` applies it: red around u' > `threshold`, blue
around u' < -`threshold` [m/s], by default 2 standard deviations of u'. The long horizontal
axis is cut to the height of the grid. Returns the `Figure`; needs a field in `am.wf`.
"""
function MakieControlPlots.plot(am::AtmosphericModel; threshold=nothing)
    field = cropped_field(am)
    threshold = something(threshold, default_threshold(field.turbulence))
    fig = Figure(; size=(900, 650))
    draw_field!(fig[1, 1], am, field, Observable(threshold))
    return fig
end

function AtmosphericModels.plot_interactive(am::AtmosphericModel)
    field = cropped_field(am)
    limit = maximum(abs, field.turbulence)
    fig = Figure(; size=(900, 700))
    slider = SliderGrid(fig[2, 1], (label="threshold", format="{:.2f} m/s",
                        range=range(limit / SLIDER_STEPS, limit; length=SLIDER_STEPS),
                        startvalue=default_threshold(field.turbulence)))
    draw_field!(fig[1, 1], am, field, only(slider.sliders).value)
    return fig
end

"""
Axes and scaled along-wind turbulence of the wind field of `am`, with the long horizontal
axis cut to the height of the grid.
"""
function cropped_field(am::AtmosphericModel)
    wf = am.wf
    wf === nothing && throw(ArgumentError(
        "no wind field: AtmosphericModel(set) only loads one when set.use_turbulence > 0"))
    x, y, z = grid_axes(am)
    height = last(z) - first(z)
    nx = findlast(<=(first(x) + height), x)
    ny = findlast(<=(first(y) + height), y)
    turbulence = wf.u[1:nx, 1:ny, :] .* turbulence_scale(am, wf)
    return (; x=x[1:nx], y=y[1:ny], z, turbulence)
end

"`THRESHOLD_SIGMAS` standard deviations of `turbulence` [m/s]."
default_threshold(turbulence) =
    THRESHOLD_SIGMAS * sqrt(sum(abs2, turbulence) / length(turbulence))

"Draw the isosurfaces of `field` at ±`threshold` [m/s] into the grid position `position`."
function draw_field!(position, am, field, threshold::Observable)
    (; x, y, z, turbulence) = field
    title = lift(threshold) do value
        "Turbulent wind field, v_wind_gnd = $(am.wf.v_wind_gnd) m/s, " *
        "threshold ± $(round(value; digits=2)) m/s"
    end
    ax = Axis3(position; aspect=:data, azimuth=-0.3π, elevation=0.2π,
               xlabel="x [m]", ylabel="y [m]", zlabel="z [m]", title)
    for (sign, color) in ((1, :red), (-1, :blue))
        contour!(ax, extrema(x), extrema(y), extrema(z), turbulence;
                 levels=lift(value -> [sign * value], threshold), colormap=[color, color],
                 alpha=ISO_ALPHA)
    end
    Legend(position, [PolyElement(; color=:red), PolyElement(; color=:blue)],
           ["gust: u' > threshold", "lull: u' < -threshold"];
           tellheight=false, tellwidth=false, halign=:right, valign=:top)
    return ax
end

end
