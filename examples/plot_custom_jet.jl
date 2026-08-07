using Pkg
if ! ("GLMakie" ∈ keys(Pkg.project().dependencies))
    Pkg.activate("examples")
end
using GLMakie, AtmosphericModels, KiteUtils

set_data_path("data")
set = load_settings("system.yaml"; relax=true)

# synthetic low-level jet: background power law plus a Gaussian bump peaking at 200 m
c_bg, a_bg, U_J, z_c, sigma = 3.0, 0.15, 4.0, 200.0, 40.0
set.heights = [10.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 600.0]
set.speeds  = [c_bg * z^a_bg + U_J * exp(-(z - z_c)^2 / (2 * sigma^2)) for z in set.heights]

heights = range(minimum(set.heights), maximum(set.heights), length=200)
v_jet = [custom_jet(set.heights, set.speeds, h) for h in heights]

fig = Figure()
ax = Axis(fig[1, 1], xlabel="height [m]", ylabel="v_wind [m/s]",
          title="Custom jet wind profile fit")
l1 = lines!(ax, heights, v_jet; linewidth=3)
sc = scatter!(ax, set.heights, set.speeds; color=:black, markersize=12)
Legend(fig[1, 2], [l1, sc], ["custom_jet", "data"])
display(fig)
