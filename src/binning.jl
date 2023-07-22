"""
Bins the data `x` and `y` using a linear binning algorithm. Returns a tuple `(g, Y, c)`

$(SIGNATURES)

## Output
- `g` - gridpoints for the x data
- `Y` - weighted average of y within each bin
- `c` - weight attached to each gridpoint

## Examples
```julia-repl
julia> x = 2π * rand(1000);

julia> y = sin.(x) + randn(size(x))/10;

julia> g, Y, c = linear_binning(x, y; nbins=100);

julia> g
100-element Vector{Float64}:
 0.007731512161049946
 0.07102439377241994
 0.13431727538378993
 ⋮
 6.210433910075309
 6.27372679168668

julia> Y
100-element Vector{Float64}:
  0.05587857685614449
  0.0670877082369906
  0.0686182710591602
  ⋮
 -0.06089998407490756
  0.007519603457183313

julia> c
100-element Vector{Float64}:
  4.7326080715176095
 14.163103524831035
  8.148664889642092
  ⋮
 11.265064011130093
  3.621345671802331
```
"""
function linear_binning(x, y; nbins=floor(Int, length(x)/100))
    g = range(minimum(x), maximum(x), length=nbins)
    grid = (g = g, c = zero(g), d = zero(g))
    linear_binning!(grid, x, y)
end

function linear_binning!(grid, x, y)
    @unpack g, c, d = grid
    for i in eachindex(x, y)
        L = (x[i] - minimum(g))/step(g) + 1 # transformation matching gridpoints to indices
        𝓁 = floor(Int, L) # index of left gridpoint
        w = 1 - (L - 𝓁) # remainder, used for weighting
        c[𝓁] += w
        d[𝓁] += w * y[i]
        d[𝓁+1] += (1-w) * y[i]
        c[𝓁+1] += 1-w
    end
    return grid
end

# TODO: combine different binning methods into one function
function simple_binning(x, y; nbins=floor(Int, length(x)/100))
    g = range(minimum(x), maximum(x), length=nbins)
    grid = (g = g, c = zero(g), d = zero(g))
    linear_binning!(grid, x, y)
end

function simple_binning!(grid, x, y)
    @unpack g, c, d = grid
    for i in eachindex(x, y)
        L = (x[i] - minimum(g))/step(g) + 1
        𝓁 = round(Int, L)
        c[𝓁] += 1
        d[𝓁] += y[i]
    end
    return grid
end
