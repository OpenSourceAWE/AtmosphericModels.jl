using AtmosphericModels

x,y,z,u,v,w,param = load_windfield(5.324)

println("Windfield dimensions: $(size(x)), $(size(y)), $(size(z))")