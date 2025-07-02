using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots, AtmosphericModels, KiteUtils

function meshgrid(x, y)
    X = [i for i in x, _ in y]
    Y = [j for _ in x, j in y]
    return (X, Y)
end

set_data_path("data") 
set = load_settings("system.yaml")
am = AtmosphericModel(set)

x,y,z,u,v,w,param = AtmosphericModels.load_windfield(am, 5.324)

X = -10:1:9
Y = -10:1:9
U, V = meshgrid(X, Y)

fig, ax = plt.subplots()
ax.quiver(X, Y, U, V)
plt.show()
