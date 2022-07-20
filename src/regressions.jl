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
    g::Vector{T}
    "Binned y data"
    Y::Vector{T}
    "Bin weights"
    c::Vector{T}

    # Pre-allocated arrays for regressions
    "Kernel weight vector"
    w::Vector{T}
    "x data centered at x₀"
    x̂::AbstractVector{T}
    "Weighting matrix"
    W::Diagonal{T, Vector{T}}
    "Polynomial basis design matrix, centered at x₀"
    X::Matrix{T}
    "`W * X`"
    WX::Matrix{T}
    "`WX' * X`"
    XWX
    "`WX' * Y`"
    XWY::Vector{T}
    "Results vector for `XWX\\XWY`"
    Σ
    ΣWX::Matrix{T}
    XWΣWX::Matrix{T}
    XWΣWXS::Matrix{T}
    V̂::Matrix{T}
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
    N = degree+1
    w = ones(T, length(g))
    W = Diagonal(w)
    X = _polybasis(g, median(g), degree)
    x̂ = view(X, :, 2)
    WX = W*X
    XWX = MMatrix{N, N, T}(WX'X)
    XWY = WX'Y

    Σ = UniformScaling(var(g))
    ΣWX = Σ*WX
    XWΣWX = WX'ΣWX
    Sₙ⁻¹ = inv(XWX)
    XWΣWXS = XWΣWX*Sₙ⁻¹
    V̂ = Sₙ⁻¹*XWΣWXS

    return LPModel{T, degree}(x, y, g, Y, c, w, x̂, W, X, WX, XWX, XWY, Σ, ΣWX, XWΣWX, XWΣWXS, V̂)
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

function _update_weights!(w, x̂, c, h; kernel=Val(:Epanechnikov))
    copyto!(w, c)
    @turbo for i in eachindex(w, x̂)
        w[i] *= Kₕ(kernel, x̂[i], h)
    end
    return w
end

function _lpreg!(::Val{N}, g, Y, c, w, x̂, W, X, WX, XWX, XWY, x₀, h; kernel=Val(:Epanechnikov)) where {N}
    _polybasis!(X, g, x₀)
    _update_weights!(w, x̂, c, h; kernel)
    mul!(WX, W, X)
    mul!(XWX, WX', X)
    mul!(XWY, WX', Y)
    # ldiv!(β, lu!(XWX), XWY)
    # return SVector{N+1, eltype(β)}(β)
    return lu(XWX)\XWY
end

function _lpreg!(𝐌::LPModel{T, N}, x₀, h; kernel=Val(:Epanechnikov)) where {T, N}
    @unpack g, Y, c, w, x̂, W, X, WX, XWX, XWY = 𝐌
    return _lpreg!(Val(N), g, Y, c, w, x̂, W, X, WX, XWX, XWY, x₀, h; kernel)
end

function _lpvcov!(ΣWX, XWΣWX, XWΣWXS, V̂, WX, XWX, Σ)
    Sₙ⁻¹ = inv(XWX)
    mul!(ΣWX, Σ, WX)
    mul!(XWΣWX, WX', ΣWX)
    rmul!(XWΣWX, 1/size(ΣWX, 1))
    mul!(XWΣWXS, XWΣWX, Sₙ⁻¹)
    mul!(V̂, Sₙ⁻¹, XWΣWXS)
    return V̂
end

function _lpvcov!(𝐌::LPModel)
    @unpack ΣWX, XWΣWX, XWΣWXS, V̂, WX, XWX, Σ = 𝐌
    return _lpvcov!(ΣWX, XWΣWX, XWΣWXS, V̂, WX, XWX, Σ)
end

function _lpvcov(::Val{N}, x̂, WX, XWX) where {N}
    σ² = var(x̂)
    Sₙ⁻¹ = inv(XWX)
    V̂ = Sₙ⁻¹ * WX' * σ² * WX * Sₙ⁻¹ / length(x̂)
    return SMatrix{N+1, N+1, eltype(V̂)}(V̂)
end

function _lpvcov(𝐌::LPModel{T, N}) where {T, N}
    @unpack x̂, WX, XWX = 𝐌
    return _lpvcov(Val(N), x̂, WX, XWX)
end

"""
Estimate the local polynomial regression model `𝐌` at the points in `v`.

$(TYPEDSIGNATURES)

## Arguments
- `𝐌::LPModel`
- `v::AbstractVector`

## Keyword Arguments
- `kernel=Val(:Epanechnikov)` - kernel function
- `h=plugin_bandwidth(x, y, size(𝐌.X, 2)-1, size(𝐌.X, 2); kernel)` - bandwidth
- `se::Bool=false` - flag for whether standard errors should be computed and returned
"""
function lpreg!(
    𝐌::LPModel{T, N},
    v::AbstractVector;
    kernel=:Epanechnikov,
    h=plugin_bandwidth(𝐌.x, 𝐌.y; ν=max(size(𝐌.X, 2)-2, 0), p=size(𝐌.X, 2)-1, kernel),
    se::Bool=false,
) where {T, N}
    # Get initial values
    β̂ = _lpreg!(𝐌, first(v), h; kernel=Val(kernel))
    se && (V̂ = SMatrix{N+1, N+1, T}(_lpvcov!(𝐌)))

    # Construct vectors for results
    𝛃 = Vector{typeof(β̂)}(undef, length(v)); 𝛃[1] = β̂
    se && (𝐕 = Vector{typeof(V̂)}(undef, length(v)); 𝐕[1] = V̂)

    # Populate vectors for remaining regressions
    for i in Base.Iterators.drop(eachindex(v), 1)
        @inbounds 𝛃[i] = _lpreg!(𝐌, v[i], h; kernel=Val(kernel))
        se && (@inbounds 𝐕[i] = SMatrix{N+1, N+1, T}(_lpvcov!(𝐌)))
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
- `kernel=Val(:Epanechnikov)` - kernel function
- `h=plugin_bandwidth(x, y, size(𝐌.X, 2)-1, size(𝐌.X, 2); kernel)` - bandwidth
- `se::Bool=false` - flag for whether standard errors should be computed and returned
"""
function lpreg(
    x::AbstractVector,
    y::AbstractVector,
    v::AbstractVector;
    degree::Int=1,
    nbins::Int=floor(Int, length(x)/100),
    kernel=:Epanechnikov,
    h=plugin_bandwidth(x, y; ν=degree-1, p=degree, kernel),
    se=false,
)
    𝐌 = LPModel(x, y, degree; nbins)
    return lpreg!(𝐌, v; kernel, h, se)
end
