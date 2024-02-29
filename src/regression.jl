struct RegressionData{T, N, P}
    y
    x
    c
    X
    x̂
    w
    W
    WX
    XWX
    XWY
end

function show(io::IO, 𝐑::RegressionData{T, N, P}) where {T, N, P}
    @unpack y, X = 𝐑
    println(io, typeof(𝐑))
    println(io, "  Variables: $N")
    println(io, "     Degree: $P")
    print(  io, "       Size: $(size(X))")
end

function RegressionData(y::Vector{T}, x::Vector{SVector{N, T}}, degree) where {N, T <: Real}
    c = ones(T, length(y))
    P = degree
    J = length(allmultiexponents(N, P))
    X = ones(T, length(x), J)
    xcols = N > 1 ? (2:1+N) : 2
    x̂ = [view(X, i, xcols) for i in axes(X, 1)]
    w = copy(vec(c))
    W = Diagonal(w)
    WX = W * X
    XWX = MMatrix{J, J}(WX'X)
    XWY = WX'y
    return RegressionData{T, N, P}(y, x, c, X, x̂, w, W, WX, XWX, XWY)
end

function RegressionData(y::Vector{T}, X::Array{T, N}, degree) where {T, N}
    x = map(SVector{size(X, 2)}, eachrow(X))
    return RegressionData(y, x, degree)
end

# function RegressionData(y::Vector{T}, X::Array{S, N}, degree) where {T <: Real, S <: Real, N}
#     R = promote_type(T, S)
#     return RegressionData(convert(Vector{R}, y), convert(Array{R, N}, X), degree)
# end

function RegressionData(grid::GridData{T, N, R}, degree) where {T, N, R}
    @unpack g, c, d = grid
    x = vec(map(SVector{N}, Base.splat(Iterators.product)(g)))
    y = copy(vec(d))
    for i in eachindex(y)
        if c[i] > 0
            y[i] /= c[i]
        end
    end
    P = degree
    J = length(allmultiexponents(N, P))
    X = ones(T, length(x), J)
    xcols = N > 1 ? (2:1+N) : 2
    x̂ = [view(X, i, xcols) for i in axes(X, 1)]
    w = copy(vec(c))
    W = Diagonal(w)
    WX = W * X
    XWX = MMatrix{J, J}(WX'X)
    XWY = WX'y
    return RegressionData{T, N, P}(y, x, vec(c), X, x̂, w, W, WX, XWX, XWY)
end

function _polybasis! end

for N in 1:5, P in 1:5
    ex = Expr(:block)
    if N > 0
        for j in 1:N
            push!(ex.args, :(X[i, $(1+j)] = xᵢ[$j] - x₀[$j]))
        end
    else
        push!(ex.args, :(X[i, 2] = xᵢ - x₀))
    end

    if P > 1
        for p in allmultiexponents(N, P) |> polypredecessors
            (i, j, k) = p
            push!(ex.args, :(X[i, $i] = X[i, $j] * X[i, $k]))
        end
    end

    @eval function _polybasis!(R::RegressionData{T, $N, $P}, X, x, x₀) where T
        for i in axes(X, 1)
            xᵢ = x[i]
            $ex
        end
        return X
    end
end

function _polybasis!(R::RegressionData, x₀)
    @unpack X, x = R
    _polybasis!(R, X, x, x₀)
end

function _update_weights!(w, x̂, c, h; kernel=Val(:Epanechnikov))
    copyto!(w, c)
    for i in eachindex(w, x̂)
        w[i] *= Kₕ(kernel, x̂[i], h)
    end
    return w
end

function locate!(R::RegressionData, x₀, h; kernel=Val(:Epanechnikov))
    @unpack X, w, x̂, c = R
    _polybasis!(R, x₀)
    _update_weights!(w, x̂, c, h; kernel)
    return R
end

function _lpreg!(𝐑::RegressionData, y, c, w, x̂, x, W, X, WX, XWX, XWY, x₀, h; kernel=Val(:Epanechnikov))
    _polybasis!(𝐑, X, x, x₀)
    _update_weights!(w, x̂, c, h; kernel)
    if sum(w) < 1e-12
        @debug "Empty neighborhood" x₀
        return nothing
    end
    mul!(WX, W, X) # including causes total 4 allocations
    mul!(XWX, WX', X) # including causes total 4 allocations
    mul!(XWY, WX', y)
    if det(XWX) < 1e-12
        @debug "Nearly singular matrix ($x₀)" XWX
        return pinv(XWX) * XWY
    else
        return lu(XWX)\XWY
    end
end

function _lpreg!(𝐑::RegressionData, x₀, h; kernel=Val(:Epanechnikov))
    @unpack y, c, w, x̂, x, W, X, WX, XWX, XWY = 𝐑
    return _lpreg!(𝐑, y, c, w, x̂, x, W, X, WX, XWX, XWY, x₀, h; kernel)
end

"""
$(SIGNATURES)

Local polynomial regression, modifying a `LPGridModel`.
"""
function lpreg!(
    𝐑::RegressionData, v;
    kernel=:Epanechnikov,
    h=plugin_bandwidth(𝐑.x, 𝐑.y; ν=max(size(𝐑.X, 2)-2, 0), p=size(𝐑.X, 2)-1, kernel)
)
    return map(x₀ -> _lpreg!(𝐑, x₀, h; kernel=Val(kernel)), v)
end
