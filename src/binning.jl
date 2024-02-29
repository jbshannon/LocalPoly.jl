"""
    togridindex(x, gmin, δ) = 1 + (x - gmin)/δ

Transform the value `x` to an index into a grid with minimum `gmin` and step length `δ`.

## Examples

```julia-repl
julia> g = 0:0.1:1
0.0:0.1:1.0

julia> gmin = minimum(g)
0.0

julia> δ = step(g)
0.1

julia> i = togridindex(0.5, gmin, δ)
6.0

julia> collect(g)[Int(i)] == 0.5
true

julia> fractional_index = togridindex(2^(-1/2), gmin, δ)
8.071067811865476
```
"""
togridindex(x, gmin, δ) = 1 + (x - gmin)/δ

"""
    togridindex(x, g::AbstractRange)

Transform the value `x` to an index into the range `g`. This method is ~3 times slower than the direct method, so repeated calls should pre-compute the range minimum and step.

## Examples

```julia-repl
julia> g = 0:0.1:1
0.0:0.1:1.0

julia> i = togridindex(0.5, g)
6.0

julia> collect(g)[Int(i)] == 0.5
true

julia> fractional_index = togridindex(2^(-1/2), gmin, δ)
8.071067811865476
```
"""
togridindex(x, g) = togridindex(x, minimum(g), step(g)) # much slower

struct GridData{T <: Real, N, R <: AbstractRange{T}}
    "Tuple of grid steps in each dimension"
    g::NTuple{N, R}

    "Weights at grid points"
    c::Array{T, N}

    "Weighted sum of data at grid points"
    d::Array{T, N}
end

function GridData(g::NTuple{N, R}) where {N, R}
    T = eltype(first(g))
    M = length.(g)
    return GridData(g, zeros(T, M...), zeros(T, M...))
end

function show(io::IO, G::GridData{T, N, R}) where {T, N, R}
    @unpack g = G
    println(io, typeof(G))
    println(io, "  Minimum: $(minimum.(g))")
    println(io, "  Maximum: $(maximum.(g))")
    println(io, "     Step: $(step.(g))")
    println(io, "   Length: $(length.(g))")
end

function show(io::IO, G::GridData{T, 1, R}) where {T, R}
    g = G.g[1]
    println(io, typeof(G))
    println(io, "  Minimum: $(minimum(g))")
    println(io, "  Maximum: $(maximum(g))")
    println(io, "     Step: $(step(g))")
    println(io, "   Length: $(length(g))")
end

# TODO: make this more sophisticated
function guessbins(X)
    N = floor(Int, size(X, 1)/10)
    return ntuple(i -> max(2, N), size(X, 2))
end

"""
$(SIGNATURES)

Linear binning algorithm for fast computation.
"""
function linear_binning(
    X::Array{T, S}, y::Vector{T};
    nbins=guessbins(X),
) where {T, S}
    N = size(X, 2)
    g = ntuple(Val(N)) do j
        Xj = view(X, :, j)
        n = nbins[j]
        range(minimum(Xj), maximum(Xj), length=n)
    end
    M = length.(g)
    grid = GridData(g, zeros(T, M), zeros(T, M))
    linear_binning!(grid, X, y)
end

function linear_binning(
    X::Vector{SVector{N, T}}, y::Vector{T};
    nbins=guessbins(X),
) where {N, T}
    g = ntuple(Val(N)) do j
        Xj = getindex.(X, j)
        n = nbins[j]
        range(minimum(Xj), maximum(Xj), length=n)
    end
    M = length.(g)
    grid = GridData(g, zeros(T, M), zeros(T, M))
    linear_binning!(grid, X, y)
end

# specialized code for 1 dimension
function linear_binning!(grid::GridData{T, 1, R}, x::Array{T, 1}, y) where {T <: Real, R}
    @unpack g, c, d = grid
    gmin = minimum(g[1])
    δ = step(g[1])
    N = last(eachindex(c))
    for i in eachindex(x, y)
        L = togridindex(x[i], gmin, δ) # transformation matching gridpoints to indices
        𝓁 = floor(Int, L) # index of left gridpoint
        w = 1 - (L - 𝓁) # remainder, used for weighting
        c[𝓁] += w
        d[𝓁] += w * y[i]
        𝓁R = min(𝓁+1, N)
        c[𝓁R] += 1-w
        d[𝓁R] += (1-w) * y[i]
    end
    return grid
end

# specialized code for 2 dimensions
function linear_binning!(
    grid::GridData{T, 2, R},
    X::Matrix{T},
    y,
) where {R, T}
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

function linear_binning!(
    grid::GridData{T, 2, R},
    X::Vector{SVector{2, T}},
    y,
) where {R, T}
    @unpack g, c, d = grid
    Ilast = last(CartesianIndices(c))
    for i in eachindex(y)
        L1 = togridindex(X[i][1], g[1])
        𝓁1 = floor(Int, L1)
        r1 = 1 - (L1 - 𝓁1)

        L2 = togridindex(X[i][2], g[2])
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

# generic code for N dimensions
"""
$(SIGNATURES)

Linear binning algorithm for fast computation.
"""
function linear_binning!(
    grid::GridData{T, N, R},
    X::Matrix{T},
    y,
) where {T, N, R}
    @unpack g = grid
    power_combinations = Base.splat(Iterators.product)(0:1 for _ in 1:N)
    ts = vec(collect(power_combinations))
    gmin = minimum.(g)
    δ = step.(g)
    for i in eachindex(y)
        _linear_binning!(grid, view(X, i, :), y[i], gmin, δ, ts)
    end
    return grid
end

function _linear_binning!(grid::GridData{T, N, R}, X, y, gmin, δ, ts) where {T, N, R}
    @unpack g, c, d = grid
    L = ntuple(j -> togridindex(X[j], gmin[j], δ[j]), N)
    𝓁 = ntuple(j -> floor(Int, L[j]), N)
    r = ntuple(j -> 1 - (L[j] - 𝓁[j]), N)
    w = ntuple(j -> prod(k -> ts[j][k] == 0 ? r[k] : 1-r[k], 1:N), 2^N)

    for j in 1:2^N
        I = CartesianIndex(𝓁) + CartesianIndex(ts[j])
        c[I] += w[j]
        d[I] += w[j] * y
    end
end

# TODO: combine different binning methods into one function
function simple_binning(x, y; nbins=guessbins(X))
    g = range(minimum(x), maximum(x), length=nbins)
    grid = (g = g, c = zero(g), d = zero(g))
    simple_binning!(grid, x, y)
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

"""
$(SIGNATURES)

Access nodes of gridded data.
"""
function gridnodes(g::NTuple{1, R}) where {R}
    return collect(first(g))
end

function gridnodes(g::NTuple{N, R}) where {N, R}
    return map(SVector{N}, Iterators.product(g...))
end

gridnodes(grid::GridData) = gridnodes(grid.g)
