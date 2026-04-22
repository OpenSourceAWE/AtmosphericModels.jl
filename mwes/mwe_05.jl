using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    Pkg.activate(joinpath(@__DIR__, "..", "examples"))
end
using ControlPlots, AtmosphericModels, KiteUtils

function meshgrid(x, y)
    X = [i for i in x, _ in y]
    Y = [j for _ in x, j in y]
    return (X, Y)
end

set_data_path("data") 
set = load_settings("system.yaml"; relax=true)
am = AtmosphericModel(set)

x,y,z,u,v,w,param = AtmosphericModels.load_windfield(am, 5.324)

x1 = x[:, 1, :][1:51,:]
y1 = y[1, :, :]
u1 = u[:, 1, :][1:51,:]
v1 = v[1, :, :]

X = -10:1:9
Y = -10:1:9
U, V = meshgrid(X, Y)

fig, ax = plt.subplots()
# ax.quiver(X, Y, U, V)
ax.quiver(x1, y1, u1, v1)
plt.show()
