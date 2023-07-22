"""
$(TYPEDEF)

---

$(TYPEDFIELDS)
"""
struct LPGridModel{T <: Real}
    # Raw data
    "Raw x data"
    x::Vector{T}
    "Raw y data"
    y::Vector{T}

    # Binned data
    "Binned x data"
    # g::Vector{T}
    g
    "Binned y data"
    Y::Vector{T}
    "Bin weights"
    c::Vector{T}

    # Pre-allocated arrays for regressions
    YΣ::Vector{T}
    WX::Matrix{T}
    S::Vector{T}
    Sₙ::AbstractMatrix{T}
    XWX
    XWY::Vector{T}
end

function LPGridModel(
    x::Vector{T},
    y::Vector{T},
    degree::Int;
    nbins::Int=floor(Int, length(x)/100),
    kernel=:Epanechnikov,
    h=plugin_bandwidth(x, y; ν=degree-1, p=degree, kernel=:Epanechnikov),
) where T <: Real

    grid = linear_binning(x, y; nbins)
    @unpack g, c, d = grid
    YΣ = d
    Δg = step(g)

    # Pre-allocate matrices
    M = 2nbins - 1
    w = zeros(T, M)
    WX = ones(T, M, 1 + 2degree)
    S = zeros(T, size(WX, 2))
    Sₙ = view(S, [i+j-1 for i in 1:1+degree, j in 1:1+degree])
    XWX = MMatrix{1+degree, 1+degree}(Sₙ)
    XWY = zeros(T, 1+degree)

    # Precompute weights and polynomial bases
    _precompute_weights!(w, Δg, h, Val(kernel))
    _precompute_polybasis!(WX, Δg)
    @turbo for j in axes(WX, 2), i in axes(WX, 1)
        WX[i, j] *= w[i]
    end

    return LPGridModel{T}(x, y, g, Y, c, YΣ, WX, S, Sₙ, XWX, XWY)
end

function _precompute_weights!(w, Δg, h, K)
    nbins = (1+length(w))÷2
    shifts = 0:nbins-1
    for n in shifts
        v = Kₕ(K, n*Δg, h)
        w[nbins+n] = v
        w[nbins-n] = v
    end
    return w
end

function _precompute_polybasis!(X, Δg)
    nbins = (1 + size(X, 1)) ÷ 2
    J = Iterators.drop(axes(X, 2), 1)

    # Compute powers of Δx
    for n in 0:nbins-1
        i = nbins+n
        Δx = n*Δg
        xᵖ = 1
        for j in J
            xᵖ *= Δx
            X[i, j] = xᵖ
        end
    end

    # Flip signs for -Δx
    for j in J, i in 1:nbins-1
        sgn = isodd(j-1) ? -1 : 1
        X[i, j] = sgn * X[end+1-i, j]
    end

    return X
end

function _shifted_matmul!(Z, gridindex, p, X, y; idxbw=length(y))
    nbins = length(y)
    fill!(Z, zero(eltype(Z)))
    I = max(gridindex-idxbw, 1):min(gridindex+idxbw, nbins)
    for j in 1:1+p, i in I
        Z[j] += X[i+nbins-gridindex, j] * y[i]
    end
    return Z
end

function _lpreg!(𝐌::LPGridModel, gridindex::Int; idxbw=length(𝐌.c))
    @unpack x, WX, c, YΣ, S, Sₙ, XWX, XWY = 𝐌
    p = (size(WX, 2)-1)÷2
    _shifted_matmul!(  S, gridindex, 2p, WX,  c; idxbw)
    _shifted_matmul!(XWY, gridindex,  p, WX, YΣ; idxbw)
    copyto!(XWX, Sₙ)
    det(XWX) < 1e-8 && @warn("Nearly singular matrix", XWX, gridindex)
    return lu(XWX)\XWY
end

function _find_idxbw(𝐌)
    x = view(𝐌.WX, : , 1)
    return (1+length(x))÷2 - findfirst(>(0), x)
end

function lpreg!(𝐌::LPGridModel; idxbw=_find_idxbw(𝐌))
    return [_lpreg!(𝐌, i; idxbw) for i in eachindex(𝐌.g)]
end

function lpreg(
    x::AbstractVector,
    y::AbstractVector;
    degree::Int=1,
    nbins::Int=floor(Int, length(x)/10),
    kernel=:Epanechnikov,
    h=plugin_bandwidth(x, y; ν=degree-1, p=degree, kernel),
    se=false,
)
    𝐌 = LPGridModel(x, y, degree; kernel, h, nbins)
    return lpreg!(𝐌)
end
