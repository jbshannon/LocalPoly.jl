kernelfunc(::Val{:Uniform}, u) = IfElse.ifelse(abs(u) <= 1, 0.5, 0.0)
kernelfunc(::Val{:Triangular}, u) = IfElse.ifelse(abs(u) <= 1, 1-abs(u), 0.0)
kernelfunc(::Val{:Epanechnikov}, u) = IfElse.ifelse(abs(u) <= 1, 3*(1-u^2)/4, 0.0)
kernelfunc(::Val{:Quartic}, u) = IfElse.ifelse(abs(u) <= 1, 15*((1-u^2)^2)/16, 0.0)
kernelfunc(::Val{:Triweight}, u) = IfElse.ifelse(abs(u) <= 1, 35*((1-u^2)^3)/32, 0.0)
kernelfunc(::Val{:Tricube}, u) = IfElse.ifelse(abs(u) <= 1, 70*((1-abs(u)^3)^3)/81, 0.0)
kernelfunc(::Val{:Cosine}, u) = IfElse.ifelse(abs(u) <= 1, (π/4)*cos((π/2)*u), 0.0)
kernelfunc(::Val{:Gaussian}, u) = pdf(Normal(), u)
kernelfunc(::Val{:Logistic}, u) = 1/(exp(u) + 2 + exp(-u))
kernelfunc(::Val{:Sigmoid}, u) = (2/π)/(exp(u) + exp(-u))
kernelfunc(::Val{:Silverman}, u) = (ū = abs(u)/√2; 0.5*exp(-ū)*sin(ū+π/4))

Kₕ(K, u::T, h::S) where {T <: Real, S <: Real} = kernelfunc(K, u/h)/h

## Multi-dimensional kernel functions
Kₕ(K, u, h) = prod(zip(u, h)) do (uu, hh)
    Kₕ(K, uu, hh)
end

# TODO: support direct computation of these constants by integrating equivalent kernels
"`Dict` containing the constant ``C_{\\nu , p}(K)`` used for the plugin bandwidth"
const 𝐶 = Dict(
    (0, 1, :Gaussian) => 0.776,
    (0, 3, :Gaussian) => 1.161,
    (1, 2, :Gaussian) => 0.884,
    (2, 3, :Gaussian) => 1.006,

    (0, 1, :Uniform) => 1.351,
    (0, 3, :Uniform) => 2.813,
    (1, 2, :Uniform) => 1.963,
    (2, 3, :Uniform) => 2.604,

    (0, 1, :Epanechnikov) => 1.719,
    (0, 3, :Epanechnikov) => 3.243,
    (1, 2, :Epanechnikov) => 2.275,
    (2, 3, :Epanechnikov) => 2.893,
    (1, 4, :Epanechnikov) => 2.8, # this is just a guess for now

    (0, 1, :Biweight) => 2.036,
    (0, 3, :Biweight) => 3.633,
    (1, 2, :Biweight) => 2.586,
    (2, 3, :Biweight) => 3.208,

    (0, 1, :Triweight) => 2.312,
    (0, 3, :Triweight) => 3.987,
    (1, 2, :Triweight) => 2.869,
    (2, 3, :Triweight) => 3.503,
)

polyderiv(β, k) = [β[n+1]*(factorial(n)÷factorial(n-k)) for n in k:length(β)-1]

"""
$(TYPEDSIGNATURES)

Estimate the rule-of-thumb plugin bandwidth.
"""
function plugin_bandwidth(
    x::AbstractVector{T}, y::AbstractVector{T};
    ν::Int=0, p::Int=ν+1, kernel=:Epanechnikov, W=LinearAlgebra.I,
) where {T <: Real}
    X = _polybasis(x, 0.0, p+3)
    β̌ = (X' * W * X)\(X' * W * y)
    σ̃² = var(y - X*β̌)
    m̌⁽ᵖ⁺¹⁾ = view(X, :, 1:3)*polyderiv(β̌, p+1)
    ȟ = 𝐶[(ν, p, kernel)] * (σ̃²/sum(abs2, m̌⁽ᵖ⁺¹⁾))^(1/(2p+3))
    return ȟ
end

function plugin_bandwidth(
    grid::GridData{T, 1, R};
    ν=0, p=ν+1, kernel=:Epanechnikov,
) where {T, R}
    @unpack g, c, d = grid
    X = _polybasis(collect(first(g)), 0.0, p+3)
    W = Diagonal(c)
    y = copy(d)
    for i in eachindex(y, c)
        if c[i] > 0
            y[i] /= c[i]
        end
    end
    β̌ = (X' * W * X)\(X' * W * y)
    σ̃² = var(y - X*β̌)
    m̌⁽ᵖ⁺¹⁾ = view(X, :, 1:3)*polyderiv(β̌, p+1)
    ȟ = 𝐶[(ν, p, kernel)] * (σ̃²/sum(abs2, m̌⁽ᵖ⁺¹⁾))^(1/(2p+3))
    return ȟ
end

## Add (approximate) finite support for kernels
support(::Val{:Uniform}) = 1
support(::Val{:Triangular}) = 1
support(::Val{:Epanechnikov}) = 1
support(::Val{:Quartic}) = 1
support(::Val{:Triweight}) = 1
support(::Val{:Tricube}) = 1
support(::Val{:Cosine}) = 1
support(::Val{:Gaussian}) = 4 # based on Wand (1997)
support(::Val{:Logistic}) = 4 # double-check the following
support(::Val{:Sigmoid}) = 4
support(::Val{:Silverman}) = 4
