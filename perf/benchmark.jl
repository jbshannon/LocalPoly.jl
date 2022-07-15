using BenchmarkTools, LocalPoly
x = 2π * rand(100_000);
y = sin.(x) + randn(size(x))/10;
v = range(minimum(x), maximum(x), length=1000);
@btime lpreg($x, $y, $v);
# 178.584 ms (6506549 allocations: 115.18 MiB)
