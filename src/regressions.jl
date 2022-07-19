"""
$(TYPEDEF)

---

$(TYPEDFIELDS)
"""
struct LPModel{T <: Real, N}
    # Raw data
    "Raw x data"
    x::AbstractVector{T}
    "Raw y data"
    y::AbstractVector{T}

    # Binned data
    "Binned x data"
    g::AbstractVector{T}
    "Binned y data"
    Y::AbstractVector{T}
    "Bin weights"
    c::AbstractVector{T}

    # Pre-allocated arrays for regressions
    "Kernel weight vector"
    w::AbstractVector{T}
    "x data centered at x₀"
    x̂::AbstractVector{T}
    "Weighting matrix"
    W::Diagonal{T, Vector{T}}
    "Polynomial basis design matrix, centered at x₀"
    X::Matrix{T}
    "`W * X`"
    WX::Matrix{T}
    "`WX' * X`"
    XWX::Matrix{T}
    "`WX' * Y`"
    XWY::Vector{T}
    "Results vector for `XWX\\XWY`"
    β::Vector{T}
end

function show(io::IO, m::LPModel)
    println(io, typeof(m))
    println(io, "        Degree: $(size(m.X, 2)-1)")
    println(io, "  Observations: $(length(m.x))")
    print(  io, "          Bins: $(length(m.g))")
end

function LPModel(x::Vector{T}, y::Vector{T}, degree::Int; nbins::Int=0) where T <: Real
    if nbins > 0
        g, Y, c = linear_binning(x, y; nbins=nbins)
    else
        g, Y, c = x, y, ones(T, size(x))
    end

    # Pre-allocate matrices
    w = ones(T, length(g))
    W = Diagonal(w)
    X = _polybasis(g, median(g), degree)
    x̂ = view(X, :, 2)
    WX = W*X
    XWX = WX'X
    XWY = WX'Y
    β = lu(XWX)\XWY

    return LPModel{T, degree}(x, y, g, Y, c, w, x̂, W, X, WX, XWX, XWY, β)
end

function LPModel(x::Vector{R}, y::Vector{S}; degree::Int=1, nbins::Int=0) where {R<:Real, S<:Real}
    x, y = promote(x, y)
    return LPModel(x, y, degree; nbins)
end

function _polybasis!(X, x, x₀)
    for i in eachindex(x)
        x̂ = x[i] - x₀
        x̂ᵖ = 1
        for p in 1:size(X, 2)-1
            x̂ᵖ *= x̂
            X[i, p+1] = x̂ᵖ
        end
    end
    return X
end

function _polybasis(x, x₀, degree)
    X = ones(eltype(x), length(x), degree+1)
    _polybasis!(X, x, x₀)
    return X
end

function _update_weights!(w, x̂, c, h; kernel=:Epanechnikov)
    Kₕ(u) = KERNELS[kernel](u/h)/h
    copyto!(w, c)
    for i in eachindex(w, x̂)
        w[i] *= Kₕ(x̂[i])
    end
    return w
end

function _lpreg!(::Val{N}, g, Y, c, w, x̂, W, X, WX, XWX, XWY, β, x₀, h; kernel=:Epanechnikov) where {N}
    _polybasis!(X, g, x₀)
    _update_weights!(w, x̂, c, h; kernel)
    mul!(WX, W, X)
    mul!(XWX, WX', X)
    mul!(XWY, WX', Y)
    ldiv!(β, lu!(XWX), XWY)
    return SVector{N+1, eltype(β)}(β)
end

function _lpreg!(𝐌::LPModel{T, N}, x₀, h; kernel=:Epanechnikov) where {T, N}
    @unpack g, Y, c, w, x̂, W, X, WX, XWX, XWY, β = 𝐌
    return _lpreg!(Val(N), g, Y, c, w, x̂, W, X, WX, XWX, XWY, β, x₀, h; kernel)
end

function _lpvcov(x̂, WX, XWX)
    σ² = var(x̂)
    Sₙ⁻¹ = inv(XWX)
    V̂ = Sₙ⁻¹ * WX' * σ² * WX * Sₙ⁻¹ / length(x̂)
    S1, S2 = size(V̂)
    return SMatrix{S1, S2, eltype(V̂)}(V̂)
end

function _lpvcov(𝐌::LPModel)
    @unpack x̂, WX, XWX = 𝐌
    return _lpvcov(x̂, WX, XWX)
end

"""
Estimate the local polynomial regression model `𝐌` at the points in `v`.

$(TYPEDSIGNATURES)

## Arguments
- `𝐌::LPModel`
- `v::AbstractVector`

## Keyword Arguments
- `kernel::Symbol=:Epanechnikov` - kernel function
- `h=plugin_bandwidth(x, y, size(𝐌.X, 2)-1, size(𝐌.X, 2); kernel)` - bandwidth
- `se::Bool=false` - flag for whether standard errors should be computed and returned
"""
function lpreg!(
    𝐌::LPModel,
    v::AbstractVector;
    kernel::Symbol=:Epanechnikov,
    h=plugin_bandwidth(𝐌.x, 𝐌.y, max(size(𝐌.X, 2)-2, 0), size(𝐌.X, 2)-1; kernel),
    se::Bool=false,
)
    # Get initial values
    β̂ = _lpreg!(𝐌, first(v), h; kernel)
    se && (V̂ = _lpvcov(𝐌))

    # Construct vectors for results
    𝛃 = repeat([β̂], length(v))
    se && (𝐕 = repeat([V̂], length(v)))

    # Populate vectors for remaining regressions
    for i in Base.Iterators.drop(eachindex(v), 1)
        @inbounds 𝛃[i] = _lpreg!(𝐌, v[i], h; kernel)
        se && (@inbounds 𝐕[i] = _lpvcov(𝐌))
    end

    return se ? (𝛃, 𝐕) : 𝛃
end

"""
Estimate the local polynomial regression of `y` on `x` at the points in `v`.

$(TYPEDSIGNATURES)

## Arguments
- `x::AbstractVector`
- `y::AbstractVector`
- `v::AbstractVector`

## Keyword Arguments
- `degree::Int=1` - degree of the polynomial approximation
- `nbins::Int=floor(Int, length(x)/100)` - number of bins to use (0 for no binning)
- `kernel::Symbol=:Epanechnikov` - kernel function
- `h=plugin_bandwidth(x, y, size(𝐌.X, 2)-1, size(𝐌.X, 2); kernel)` - bandwidth
- `se::Bool=false` - flag for whether standard errors should be computed and returned
"""
function lpreg(
    x::AbstractVector,
    y::AbstractVector,
    v::AbstractVector;
    degree::Int=1,
    nbins::Int=floor(Int, length(x)/100),
    kernel::Symbol=:Epanechnikov,
    h=plugin_bandwidth(x, y, degree-1, degree; kernel),
    se=false,
)
    𝐌 = LPModel(x, y, degree; nbins)
    return lpreg!(𝐌, v; kernel, h, se)
end
