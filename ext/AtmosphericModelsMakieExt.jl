module AtmosphericModelsMakieExt

using MakieControlPlots.Makie
import MakieControlPlots
using AtmosphericModels
using AtmosphericModels: grid_axes

"""
    plot(am::AtmosphericModel)

Plot the turbulent wind field of `am` in 3D: the along-wind turbulence [m/s], scaled as
[`get_wind`](@ref) applies it, on the three visible faces of the grid. The long horizontal
axis is cut to the height of the grid. Returns the `Figure`; needs a field in `am.wf`.
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
    rel_turb = am.set.use_turbulence * rel_turbo(am, wf.v_wind_gnd)
    turbulence = wf.u[1:nx, 1:ny, :] .* rel_turb
    limit = maximum(abs, turbulence)

    fig = Figure(; size=(900, 650))
    ax = Axis3(fig[1, 1]; aspect=:data, azimuth=-0.3π, elevation=0.2π,
               xlabel="x [m]", ylabel="y [m]", zlabel="z [m]",
               title="Turbulent wind field, v_wind_gnd = $(wf.v_wind_gnd) m/s")
    style = (; colormap=:balance, colorrange=(-limit, limit), shading=NoShading)
    Y, Z = face(y, z)
    downwind = surface!(ax, fill(x[end], size(Y)), Y, Z; color=turbulence[end, :, :],
                        style...)
    X, Z = face(x, z)
    surface!(ax, X, fill(y[1], size(X)), Z; color=turbulence[:, 1, :], style...)
    X, Y = face(x, y)
    surface!(ax, X, Y, fill(z[end], size(X)); color=turbulence[:, :, end], style...)
    Colorbar(fig[1, 2], downwind; label="along-wind turbulence u' [m/s]")
    return fig
end

"Coordinate matrices of the grid face spanned by the axes `first_axis` and `second_axis`."
face(first_axis, second_axis) =
    (first_axis .* one.(second_axis)', one.(first_axis) .* second_axis')

end
