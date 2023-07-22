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
function linear_binning(x, y::Vector{T}; nbins=floor(Int, length(x)/100)) where {T}
    g = range(minimum(x), maximum(x), length=nbins)
    M = length(g)
    grid = (g = g, c = zeros(T, M), d = zeros(T, M))
    linear_binning!(grid, x, y)
end

function linear_binning(X::Matrix{T}, y; nbins=ntuple(i -> floor(Int, size(X, 1)/100), size(X, 2))) where {T}
    g = ntuple(size(X, 2)) do j
        Xj = view(X, :, j)
        range(minimum(Xj), maximum(Xj), length=nbins[j])
    end
    M = length.(g)
    grid = (g = g, c = zeros(T, M), d = zeros(T, M))
    linear_binning!(grid, x, y)
end

# generic code for N dimensions
function _linear_binning!(grid, X, y, ts, ::Val{N}) where N
    @unpack g, c, d = grid
    L = ntuple(j -> 1 + (X[j] - minimum(g[j]))/step(g[j]), Val(N))
    𝓁 = ntuple(j -> floor(Int, L[j]), Val(N))
    r = ntuple(j -> 1 - (L[j] - 𝓁[j]), Val(N))
    w = ntuple(j -> prod(k -> ts[j][k] == 0 ? r[k] : 1-r[k], 1:N), Val(2^N))

    for j in 1:2^N
        I = CartesianIndex(𝓁) + CartesianIndex(ts[j])
        c[I] += w[j]
        d[I] += w[j] * y
    end
end

# specialized code for 1 dimension
function linear_binning!(grid, x::Array{T, 1}, y) where {T <: Real}
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

# specialized code for 2 dimensions
function linear_binning!(grid, X::Array{T, 2}, y) where {T <: Real}
    @unpack g, c, d = grid
    for i in eachindex(y)
        L1 = (X[i, 1] - gmin[1])/δ[1] + 1
        𝓁1 = floor(Int, L1)
        r1 = 1 - (L1 - 𝓁1)

        L2 = (X[i, 2] - gmin[2])/δ[2] + 1
        𝓁2 = floor(Int, L2)
        r2 = 1 - (L2 - 𝓁2)

        w = r1 * r2
        c[𝓁1, 𝓁2] += w
        d[𝓁1, 𝓁2] += w * y[i]

        w = (1 - r1) * r2
        c[𝓁1+1, 𝓁2] += w
        d[𝓁1+1, 𝓁2] += w * y[i]

        w = r1 * (1 - r2)
        c[𝓁1, 𝓁2+1] += w
        d[𝓁1, 𝓁2+1] += w * y[i]

        w = (1 - r1) * (1 - r2)
        c[𝓁1+1, 𝓁2+1] += w
        d[𝓁1+1, 𝓁2+1] += w * y[i]
    end
    return c, d
end

@generated function linear_binning!(grid, X::Array{T, N}, y) where {T <: Real, N}
    quote
        power_combinations = @ncall $N Iterators.product i -> 0:1
        ts = vec(collect(power_combinations))
        _linear_binning!(grid, X, y, ts, Val(N))
    end
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
