using BenchmarkTools, LocalPoly
x = 2π * rand(100_000)
y = sin.(x) + randn(size(x))/10
v = range(minimum(x), maximum(x), length=1000)
@btime h = plugin_bandwidth($x, $y)
# 2.701 ms (14 allocations: 6.10 MiB)
@btime lpreg($x, $y, $v; h=$h)
# 178.584 ms (6506549 allocations: 115.18 MiB)
# 11.204 ms (9563 allocations: 442.02 KiB)

## Time LPModel creation
@btime 𝐌 = LPModel($x, $y, 1; nbins=1000)
# 10.051 ms (9491 allocations: 351.53 KiB)

## Time a single regression
@btime LocalPoly._lpreg!($𝐌, $v[end÷2], $h);
# 9.346 μs (6 allocations: 160 bytes)
