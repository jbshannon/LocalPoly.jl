# LocalPoly.jl

[![Build Status](https://github.com/jbshannon/LocalPoly.jl/workflows/CI/badge.svg)](https://github.com/jbshannon/LocalPoly.jl/actions?query=workflows/CI)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://jbshannon.github.io/LocalPoly.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://jbshannon.github.io/LocalPoly.jl/dev)

`LocalPoly.jl` is a Julia implementation of the local polynomial regression methods outlined in [Fan and Gijbels (1996)](https://doi.org/10.1201/9780203748725). This package is still experimental, and the API is subject to change.

## Overview

This package provides the function `lpreg`, which computes the local polynomial regression coefficients for input data of any dimension. Since the local polynomial estimator is much faster over a grid of evenly-spaced data points, performance can be improved by first projecting the data to such a grid using the `linear_binning` function.

## Examples

```julia
using LocalPoly, Random
Random.seed!(42)
x = 2π * rand(1000)
y = sin.(x) + randn(size(x))/4
grid = linear_binning(x, y; nbins=100)
β̂ = lpreg(grid)
```

The first element of the coefficient vector represents the function estimate at the point of evaluation:

```julia
julia> ŷ = first.(β̂)
100-element Vector{Float64}:
 -0.03070776997429395
  0.048352477003287916
  ⋮
 -0.04452583837750935
 -0.04543586963674676
```

The estimation points can be recovered by using the `gridnodes` method:

```julia
julia> v = gridnodes(grid)
100-element Vector{Float64}:
 0.0005470440483675072
 0.06395994969485173
 ⋮
 6.215011797403821
 6.278424703050305
```

Plotting the fitted function values against the data:

```julia
using CairoMakie
f = Figure()
ax = Axis(f[1, 1])
scatter!(ax, x, y; markersize=2, label="Data")
lines!(ax, v, sin.(v); color=Cycled(2), label="True values")
lines!(ax, v, ŷ; color=:tomato, linewidth=3, label="Fitted values")
Legend(f[2, 1], ax; orientation=:horizontal, framevisible=false)
current_figure()
```

![Fit](./docs/src/images/readme/light/fit.svg#gh-light-mode-only)
![Fit](./docs/src/images/readme/dark/fit.svg#gh-dark-mode-only)

### Standard Errors

The bias and variance of the estimated coefficients can be estimated by fitting a "pilot" 
model using a different bandwidth and higher-degree polynomial approximation (see
documentation for details):

```julia
julia> ci = confint(grid; α=0.05)
100-element Vector{Tuple{Float64, Float64}}:
 (-0.01564668938703559, 0.2698268678178526)
 (0.04315971509904605, 0.22372303301848356)
 (0.09309771887828536, 0.23823437533995848)
 ⋮
 (-0.24787454088605573, -0.024976639691117714)
 (-0.2363528259834405, 0.04772308836500469)
 (-0.3155131059523961, 0.12284688160121504)
```

Adding the confidence intervals to the plot:

```julia
band!(ax, v, first.(ci), last.(ci); color=(:tomato, 0.3))
current_figure()
```

![Fit with confidence interval](./docs/src/images/readme/light/fit_ci.svg#gh-light-mode-only)
![Fit with confidence interval](./docs/src/images/readme/dark/fit_ci.svg#gh-dark-mode-only)

### Derivative Estimation

One of the advantages of the local polynomial estimator is its ability to nonparametrically estimate the derivatives of the regression function. To estimate the $\nu$-th derivative of the function, we simply fit a model with degree of at least $\nu+1$.

```julia
ν = 1
degree = ν+1
h = plugin_bandwidth(grid; ν)
β̂1 = lpreg(grid; degree, h)
ŷ′ = [b[ν+1] for b in β̂1]
ci1 = confint(grid; ν, p=degree)

fig, ax, sc = scatter(x, y; color=Cycled(1), markersize=3, label="Data")
lines!(ax, v, cos.(v); label="True values")
lines!(ax, v, ŷ′; linewidth=3, label="Fitted values")
band!(ax, v, first.(ci1), last.(ci1); color=(:tomato, 0.3))
Legend(fig[2, 1], ax; orientation=:horizontal)
current_figure()
```

![Fit derivative curve](./docs/src/images/readme/light/fit_derivative.svg#gh-light-mode-only)
![Fit derivative curve](./docs/src/images/readme/dark/fit_derivative.svg#gh-dark-mode-only)

Fits can be improved with more sophisticated bandwidth selection techniques (TBD).

## Performance

Set the number of observations to 100,000 and $Y_i = \sin(X_i) + \varepsilon_i$ for $X_i \in [0, 2\pi]$. Evaluate the local polynomial estimator at 1,000 points.

```julia
using BenchmarkTools, LocalPoly
x = 2π * rand(100_000)
y = sin.(x) + randn(size(x))/10
@btime grid = linear_binning($x, $y; nbins=1000)
# 1.072 ms (2 allocations: 15.88 KiB)
@btime h = plugin_bandwidth($grid)
# 5.825 μs (14 allocations: 14.23 KiB)
@btime lpreg($grid; h=$h)
# 136.371 μs (392 allocations: 35.23 KiB)
```

### R

```R
library(KernSmooth)
library(microbenchmark)
x <- 2*pi*runif(100000)
y <- sin(x) + rnorm(100000)/10
v <- seq(from = 0, to = 2*pi, length.out = 1000)
h <- dpill(x, y, gridsize = 1000, range.x = c(0, 2*pi))
microbenchmark("KernSmooth" = locpoly(x, y, bandwidth = h, gridsize = 1000, range.x = c(0, 2*pi)))
```
Output:
```
Unit: milliseconds
       expr      min       lq     mean   median       uq      max neval
 KernSmooth 2.062024 2.992719 3.506988 3.205222 3.713487 12.05903   100
```

### Stata

```stata
clear all
qui set obs 100000
gen x = 2*3.14159265*runiform()
gen y = sin(x) + rnormal()/10
forval i = 1/10 {
    timer on `i'
    lpoly y x, n(1000) kernel(epan2) degree(1) nograph
    timer off `i'
}
timer list
```

Output (measured in seconds):
```
1:     14.59 /        1 =      14.5850
2:     14.45 /        1 =      14.4500
3:     14.07 /        1 =      14.0730
4:     14.31 /        1 =      14.3090
5:     14.44 /        1 =      14.4440
6:     14.31 /        1 =      14.3120
7:     14.06 /        1 =      14.0630
8:     14.22 /        1 =      14.2160
9:     14.33 /        1 =      14.3280
10:     15.00 /        1 =      14.9980
```

### MATLAB

```matlab
x = rand(100000, 1);
y = sin(x) + randn(100000, 1)/10;
v = linspace(min(x), max(x), 1000);
h = 0.087; % approximate plugin bandwidth

T = 100; % number of benchmark trials to run
tocs = zeros(T, 1);
for i = 1:numel(tocs)
    tic; lpreg(x, y, v, h); tocs(i) = toc;
end
fprintf('mean = %4.3f s\n std = %4.3f s\n', mean(tocs), std(tocs));

function betas = lpreg(x, y, v, h)
    X = [ones(size(x)) x];
    betas = v;
    for i = 1:numel(v)
        d = x - v(i);
        X(:, 2) = d;
        w = kernelfunc(d/h)/h;
        beta = inv(X' * (w .* X))*(X' * (w .* y));
        betas(i) = beta(1);
    end

    function z = kernelfunc(u)
        I = abs(u) <= 1;
        z = zeros(size(u));
        z(I) = 3*(1-u(I).^2)/4;
    end
end
```
Output:
```
mean = 2.739 s
 std = 0.130 s
```

## References

1. Fan, J., & Gijbels, I. (1996). Local Polynomial Modelling and its Applications (1st ed.). Chapman & Hall.
2. Wand, M. P. (1994). Fast Computation of Multivariate Kernel Estimators. Journal of Computational and Graphical Statistics, 3(4), 433–445. https://doi.org/10.2307/1390904
