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
    xmin, xmax = extrema(x)
    Δ = abs(xmax - xmin)/(nbins-1) # ensure maximum(x) < maximum(g)
    g = collect(range(xmin, step=Δ, length=nbins))
    Y = zeros(eltype(g), nbins)
    c = zeros(eltype(g), nbins)
    linear_binning!(g, Y, c, x, y)
end

function linear_binning!(g, Y, c, x, y)
    xmin, xmax = extrema(x)
    Δ = abs(xmax - xmin)/(length(g)-1) # ensure maximum(x) < maximum(g)

    @turbo for i in eachindex(x, y)
        L = (x[i] - xmin)/Δ + 1 # transformation matching gridpoints to indices
        𝓁 = floor(Int, L) # index of left gridpoint
        w = L - 𝓁 # remainder, used for weighting
        Y[𝓁] += w * y[i]
        Y[𝓁+1] += (1-w) * y[i]
        c[𝓁] += w
        c[𝓁+1] += 1-w
    end

    # Replace weighted sum in Y with weighted average
    for i in eachindex(Y, c)
        if c[i] > 0
            @inbounds Y[i] /= c[i]
        end
    end

    return g, Y, c
end

# TODO: combine different binning methods into one function
function simple_binning(x, y; nbins=floor(Int, length(x)/100))
    xmin, xmax = extrema(x)
    Δ = abs(xmax - xmin)/(nbins - 1) # ensure maximum(x) < maximum(g)
    g = collect(range(xmin, step=Δ, length=nbins))

    Y = zeros(eltype(g), nbins)
    c = zeros(eltype(g), nbins)
    @turbo for i in eachindex(x, y)
        L = (x[i] - xmin)/Δ + 1
        𝓁 = round(Int, L)
        Y[𝓁] += y[i]
        c[𝓁] += 1
    end

    # replace weighted sum in Y with weighted average
    for i in eachindex(Y, c)
        if c[i] > 0
            Y[i] /= c[i]
        end
    end

    return g, Y, c
end
