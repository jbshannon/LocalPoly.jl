"""
    togridindex(x, g)

Transform the value `x` to an index into the range `g`.
"""
togridindex(x, g) = 1 + (x - minimum(g))/step(g)

"""
Bins the data `x` and `y` using a linear binning algorithm.

$(SIGNATURES)

## Output
A `NamedTuple` with fields:
- `g` - `Tuple` of `StepRange`s representing the grid points
- `c` - weight attached to each gridpoint
- `d` - weighted sum of `y` at each gridpoint

## Examples
```julia-repl
julia> linear_binning([0, 0.3, 1], [1, 1, 1]*2; nbins=2)
(g = 0.0:1.0:1.0, c = [1.7, 1.3], d = [3.4, 2.6])
```
"""
function linear_binning(
    x::Vector{T}, y::Vector{S};
    nbins=max(2, floor(Int, length(x)/100))
) where {T <: Real, S <: Real}
    g = range(minimum(x), maximum(x), length=nbins)
    R = promote_type(typeof(togridindex(first(x), g)), S)
    M = length(g)
    grid = (g = g, c = zeros(R, M), d = zeros(R, M))
    linear_binning!(grid, x, y)
end

function linear_binning(
    X::Matrix{T}, y::Vector{T};
    nbins=ntuple(i -> max(2, floor(Int, size(X, 1)/100)), size(X, 2))
) where {T <: Real}
    g = ntuple(size(X, 2)) do j
        Xj = view(X, :, j)
        range(minimum(Xj), maximum(Xj), length=nbins[j])
    end
    M = length.(g)
    grid = (g = g, c = zeros(T, M), d = zeros(T, M))
    linear_binning!(grid, X, y)
end

# generic code for N dimensions
function _linear_binning!(grid, X, y, ts, ::Val{N}) where N
    @unpack g, c, d = grid
    L = ntuple(j -> togridindex(X[j], g[j]), Val(N))
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
    Ilast = last(CartesianIndices(c))
    for i in eachindex(x, y)
        L = togridindex(x[i], g) # transformation matching gridpoints to indices
        𝓁 = floor(Int, L) # index of left gridpoint
        w = 1 - (L - 𝓁) # remainder, used for weighting
        c[𝓁] += w
        d[𝓁] += w * y[i]
        I𝓁 = min(CartesianIndex(𝓁 + 1), Ilast)
        c[I𝓁] += 1-w
        d[I𝓁] += (1-w) * y[i]
    end
    return grid
end

# specialized code for 2 dimensions
function linear_binning!(
    grid::NamedTuple{(:g, :c, :d), Tuple{NTuple{2, R}, Array{T, 2}, Array{T, 2}}},
    X::Matrix{T},
    y,
) where {R, T <: Real}
    @unpack g, c, d = grid
    Ilast = last(CartesianIndices(c))
    for i in eachindex(y)
        L1 = togridindex(X[i, 1], g[1])
        𝓁1 = floor(Int, L1)
        r1 = 1 - (L1 - 𝓁1)

        L2 = togridindex(X[i, 2], g[2])
        𝓁2 = floor(Int, L2)
        r2 = 1 - (L2 - 𝓁2)

        w = r1 * r2
        c[𝓁1, 𝓁2] += w
        d[𝓁1, 𝓁2] += w * y[i]

        w = (1 - r1) * r2
        I𝓁 = min(CartesianIndex(𝓁1+1, 𝓁2), Ilast)
        c[I𝓁] += w
        d[I𝓁] += w * y[i]

        w = r1 * (1 - r2)
        I𝓁 = min(CartesianIndex(𝓁1, 𝓁2+1), Ilast)
        c[I𝓁] += w
        d[I𝓁] += w * y[i]

        w = (1 - r1) * (1 - r2)
        I𝓁 = min(CartesianIndex(𝓁1+1, 𝓁2+1), Ilast)
        c[I𝓁] += w
        d[I𝓁] += w * y[i]
    end
    return grid
end

@generated function linear_binning!(
    grid::NamedTuple{(:g, :c, :d), Tuple{NTuple{N, R}, Array{T, N}, Array{T, N}}},
    X::Matrix{T},
    y,
) where {N, R, T <: Real}
    quote
        power_combinations = @ncall $N Iterators.product i -> 0:1
        ts = vec(collect(power_combinations))
        for i in eachindex(y)
            _linear_binning!(grid, view(X, i, :), y[i], ts, Val(N))
        end
        return grid
    end
end

# TODO: combine different binning methods into one function
function simple_binning(x, y; nbins=max(2, floor(Int, length(x)/100)))
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
