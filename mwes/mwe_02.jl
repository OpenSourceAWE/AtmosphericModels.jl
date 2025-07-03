using BenchmarkTools
const A=rand(2026,51,251)
@btime maximum(A)
#   6.946 ms (0 allocations: 0 bytes) # 4.569 ms (0 allocations: 0 bytes)
# @btime treduce(max, A)
# 3.589 ms (126 allocations: 9.86 KiB)
@btime @fastmath maximum(A)
# 4.444 ms on battery