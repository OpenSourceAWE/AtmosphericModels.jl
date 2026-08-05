using Pkg
if ! ("GLMakie" ∈ keys(Pkg.project().dependencies))
    Pkg.activate(joinpath(@__DIR__, "..", "examples"))
end
using GLMakie

function meshgrid(x, y)
    X = [i for i in x, _ in y]
    Y = [j for _ in x, j in y]
    return (X, Y)
end

X = -10:1:9
Y = -10:1:9
U, V = meshgrid(X, Y)

fig = Figure()
ax = Axis(fig[1, 1])
arrows2d!(ax, vec(X), vec(Y), vec(U), vec(V))
display(GLMakie.Screen(), fig)
