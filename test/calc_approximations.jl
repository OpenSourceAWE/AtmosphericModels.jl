using Pkg
if ! ("Remez" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using AtmosphericModels, BenchmarkTools, Remez, KiteUtils

set_data_path()
set = load_settings("system.yaml"; relax=true)
am::AtmosphericModel = AtmosphericModel(set)

tofloat64(x) = [Float64(elem) for elem in x]
function wind_factor1(height)
    calc_wind_factor(am, height, Val{1})
end
function wind_factor2(height)
    calc_wind_factor(am, height, Val{2})
end
function wind_factor3(height)
    calc_wind_factor(am, height, Val{3})
end

res1=ratfn_minimax(x->wind_factor1(1/x), (1/1000,1/6), 10, 2)
res2=ratfn_minimax(x->wind_factor2(1/x), (1/1000,1/6), 10, 2)
res3=ratfn_minimax(x->wind_factor3(1/x), (1/1000,1/6), 10, 2)

println("profile_law 1:")
println(Tuple(tofloat64(res1[1])))
println(Tuple(tofloat64(res1[2])))
println("profile_law 2:")
println(Tuple(tofloat64(res2[1])))
println(Tuple(tofloat64(res2[2])))
println("profile_law 3:")
println(Tuple(tofloat64(res3[1])))
println(Tuple(tofloat64(res3[2])))
