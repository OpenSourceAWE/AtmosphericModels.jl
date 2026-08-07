using Pkg
if ! ("GLMakie" ∈ keys(Pkg.project().dependencies))
    Pkg.activate("examples")
end
using GLMakie, AtmosphericModels, KiteUtils

set_data_path("data")
set = load_settings("system.yaml"; relax=true)
set.heights = [6.0, 20.0, 60.0, 100.0, 200.0]
set.speeds  = [5.324, 6.8, 8.1, 8.9, 9.8]

heights = range(minimum(set.heights), maximum(set.heights), length=200)
v_log = [custom_log(set.heights, set.speeds, h) for h in heights]
v_exp = [custom_exp(set.heights, set.speeds, h) for h in heights]

fig = Figure()
ax = Axis(fig[1, 1], xlabel="height [m]", ylabel="v_wind [m/s]",
          title="Custom log/exp wind profile fit")
l1 = lines!(ax, heights, v_log; linewidth=3)
l2 = lines!(ax, heights, v_exp; linewidth=3)
sc = scatter!(ax, set.heights, set.speeds; color=:black, markersize=12)
Legend(fig[1, 2], [l1, l2, sc], ["custom_log", "custom_exp", "data"])
display(fig)
