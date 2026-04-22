using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    Pkg.activate(joinpath(@__DIR__, "..", "examples"))
end
using ControlPlots

function meshgrid(x, y)
    X = [i for i in x, _ in y]
    Y = [j for _ in x, j in y]
    return (X, Y)
end

X = -10:1:9
Y = -10:1:9
U, V = meshgrid(X, Y)

fig, ax = plt.subplots()
q = ax.quiver(X, Y, U, V)
# ax.quiverkey(q, X=0.3, Y=1.1, U=10, label="Quiver key, length = 10", labelpos="E")

plt.show()
