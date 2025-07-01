using AtmosphericModels

wf = load_windfield(8.0)
x, y, z = wf.x, wf.y, wf.z
println("Windfield dimensions: $(size(x)), $(size(y)), $(size(z))")