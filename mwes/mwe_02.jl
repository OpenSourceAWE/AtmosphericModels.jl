using BenchmarkTools
const A=rand(2026,51,251)
@btime maximum(A)
#   6.946 ms (0 allocations: 0 bytes)