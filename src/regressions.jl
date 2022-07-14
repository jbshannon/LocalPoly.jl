struct LPModel{T <: Real}
    # Raw data
    x::AbstractVector{T}
    y::AbstractVector{T}

    # Binned data
    g::AbstractVector{T}
    Y::AbstractVector{T}
    c::AbstractVector{T}

    # Pre-allocated arrays for regressions
    w::AbstractVector{T}
    x̂::AbstractVector{T}
    W::Diagonal{T, Vector{T}}
    X::Matrix{T}
    WX::Matrix{T}
    XWX::Matrix{T}
    XWY::Vector{T}
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
    N = length(g)
    k = degree+1
    w = zeros(T, N)
    W = Diagonal(w)
    X = ones(T, N, k)
    x̂ = view(X, :, 2)
    WX = ones(T, N, k)
    XWX = ones(T, k, k)
    XWY = ones(T, k)

    return LPModel(x, y, g, Y, c, w, x̂, W, X, WX, XWX, XWY)
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
    K = KERNELS[kernel]
    copyto!(w, c)
    for i in eachindex(w, x̂)
        w[i] *= h * K(x̂[i]/h)
    end
    return w
end

function _lpreg!(g, Y, c, w, x̂, W, X, WX, XWX, XWY, h, x₀; kernel=:Epanechnikov)
    _polybasis!(X, g, x₀)
    _update_weights!(w, x̂, c, h; kernel)
    mul!(WX, W, X)
    mul!(XWX, WX', X)
    mul!(XWY, WX', Y)
    β̂ = XWX\XWY
    return SVector{length(β̂), eltype(β̂)}(β̂)
end

function _lpreg!(𝐌::LPModel, h, x₀; kernel=:Epanechnikov)
    @unpack g, Y, c, w, x̂, W, X, WX, XWX, XWY = 𝐌
    return _lpreg!(g, Y, c, w, x̂, W, X, WX, XWX, XWY, h, x₀; kernel)
end

function _lpvcov(x̂, WX, XWX)
    σ² = var(x̂)
    Sₙ⁻¹ = inv(XWX)
    V̂ = Sₙ⁻¹ * WX' * σ² * WX * Sₙ⁻¹
    S1, S2 = size(V̂)
    return SMatrix{S1, S2, eltype(V̂)}(V̂)
end

function _lpvcov(𝐌::LPModel)
    @unpack x̂, WX, XWX = 𝐌
    return _lpvcov(x̂, WX, XWX)
end

function lpreg(
    x::AbstractVector, y::AbstractVector, v::AbstractVector;
    ν::Int=0,
    degree::Int=1,
    nbins::Int=0,
    kernel::Symbol=:Epanechnikov,
    h=plugin_bandwidth(x, y, ν, degree; kernel),
    se=false,
)
    𝐌 = LPModel(x, y; degree, nbins)

    # Get initial values
    β̂ = _lpreg!(𝐌, h, first(v); kernel)
    se && (V̂ = _lpvcov(𝐌))

    # Construct vectors for results
    𝛃 = repeat([β̂], length(v))
    se && (𝐕 = repeat([V̂], length(v)))

    # Populate vectors for remaining regressions
    for i in Base.Iterators.drop(eachindex(v), 1)
        @inbounds 𝛃[i] = _lpreg!(𝐌, h, v[i]; kernel)
        se && (@inbounds 𝐕[i] = _lpvcov(𝐌))
    end

    return se ? (𝛃, 𝐕) : 𝛃
end
