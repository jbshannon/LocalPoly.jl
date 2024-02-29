# TODO: implement convolution form

function _bias(𝐑::RegressionData{T, 1, P}, β̂, p, h; kernel=:Epanechnikov) where {T, P}
    @unpack x, X, W, WX, XWX = 𝐑
    β = ones(T, size(β̂[1]))
    βₚ = view(β, 2+p:length(β))
    Xₐ = view(X, :, 1:1+p)
    Xₚ = view(X, :, 2+p:size(X, 2))
    WXₐ = W * Xₐ
    # WXₚ = W * Xₚ
    XWXₐ = MMatrix{1+p, 1+p}(WXₐ'Xₐ)
    τ = Xₚ * βₚ
    XWτ = WXₐ'τ
    b̂ = map(eachindex(x, β̂)) do i
        locate!(𝐑, x[i], h; kernel=Val(kernel))
        copyto!(β, β̂[i])
        mul!(τ, Xₚ, βₚ)
        mul!(WXₐ, W, Xₐ)
        mul!(XWτ, WXₐ', τ)
        mul!(XWXₐ, WXₐ', Xₐ)
        return lu(XWXₐ)\XWτ
    end
    return b̂
end

function _std(𝐑::RegressionData, β̂, h; kernel=:Epanechnikov)
    @unpack y, x, X, w, W, WX, XWX = 𝐑
    XWWX = WX'WX
    ŷ = zero(y)
    map(eachindex(x, β̂)) do i
        num = zero(eltype(y))
        locate!(𝐑, x[i], h; kernel=Val(kernel))
        mul!(ŷ, X, β̂[i])
        @inbounds for j in eachindex(y, ŷ, w)
            num += (y[j] - ŷ[j])^2 * w[j]
        end

        mul!(WX, W, X)
        mul!(XWX, WX', X)
        mul!(XWWX, WX', WX)
        den = tr(W) - tr(XWWX * inv(XWX))

        return num/den
    end
end

function _var(𝐑::RegressionData, β̂, p, h; kernel=:Epanechnikov)
    @unpack x, X, W, WX, XWX = 𝐑
    σ̂² = _std(𝐑, β̂, h; kernel)
    Xₐ = view(X, :, 1:1+p)
    WXₐ = W * Xₐ
    XWXₐ = MMatrix{1+p, 1+p}(WXₐ'Xₐ)
    XWWXₐ = MMatrix{1+p, 1+p}(WXₐ'WXₐ)
    V̂ = map(eachindex(x, β̂, σ̂²)) do i
        locate!(𝐑, x[i], h; kernel=Val(kernel))
        mul!(WXₐ, W, Xₐ)
        mul!(XWXₐ, WXₐ', Xₐ)
        mul!(XWWXₐ, WXₐ', WXₐ)
        XWXinv = inv(XWXₐ)
        return XWXinv * XWWXₐ * XWXinv * σ̂²[i]
    end
    return V̂
end

function _tcrit(𝐑::RegressionData, h; α=0.05, kernel=:Epanechnikov)
    @unpack x, w = 𝐑
    crit = 1-α/2
    map(eachindex(x, w)) do i
        locate!(𝐑, x[i], h; kernel=Val(kernel))
        dₙ = floor(Int, (sum(w)^2)/sum(abs2, w))
        return quantile(TDist(dₙ), crit)
    end
end

"""
$(SIGNATURES)

Generate a confindence interval.
"""
function confint(
    grid::GridData{T, 1, R};
    ν=0, p=ν+1, a=2, α=0.05,
    kernel=:Epanechnikov,
    hpilot=plugin_bandwidth(grid; ν, p=p+a, kernel),
    h=plugin_bandwidth(grid; ν, p, kernel),
) where {T, R}
    β̂pilot = lpreg(grid; degree=p+a, h=hpilot, kernel)
    𝐑pilot = RegressionData(grid, p+a)
    b̂ = _bias(𝐑pilot, β̂pilot, p, h; kernel)
    V̂ = _var(𝐑pilot, β̂pilot, p, h; kernel)
    t = _tcrit(𝐑pilot, hpilot; α, kernel)

    β̂ = lpreg(grid; degree=p, h=h, kernel)
    𝐑 = RegressionData(grid, p)
    @unpack x, w = 𝐑
    ν! = factorial(ν)
    map(eachindex(β̂)) do i
        locate!(𝐑, x[i], h; kernel=Val(kernel))
        B̂ = sum(j -> b̂[j][ν+1] * w[j], eachindex(b̂, w))/sum(w)
        # B̂ = b̂[i][ν+1]
        center = β̂[i][ν+1] - ν! * B̂
        shift = ν! * t[i] * sqrt(V̂[i][ν+1, ν+1])
        return (center-shift, center+shift)
    end
end
