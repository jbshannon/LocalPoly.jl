module LocPolyReg

using Distributions
using StatsBase
using LinearAlgebra
import Polynomials

"Build the pseudo-Vandermonde design matrix"
function build_X₀!(X₀, x̄, 𝑝)
    # Fill matrix
    if 𝑝 >= 1
        X₀[:, 2] = x̄
    end

    if 𝑝 >= 2
        for p = 2:𝑝
            X₀[:, p+1] = x̄.^p
        end
    end

    return X₀
end

function build_X₀(x̄, 𝑝)
    X₀ = ones(length(x̄), 𝑝+1)
    build_X₀!(X₀, x̄, 𝑝)
    return X₀
end

"Fan-Gijbels plugin bandwidth estimator"
function plugin_bandwidth(x, y)
    𝐩 = Polynomials.fit(x, y, 5)
    𝐩″ = Polynomials.derivative(𝐩, 2)
    ε̂ = y .- 𝐩.(x)
    σ̃² = var(ε̂)
    h = 2.275 * (σ̃² / sum(abs2∘𝐩″, x))^(1/7)
    return h
end

function local_poly_reg(x, y, x₀, h)
    x̄ = x .- x₀
    X₀ = build_X₀(x̄, 2)
    K(u) = pdf(Epanechnikov(), u)
    W = Diagonal(K.(x̄ ./ h))
    β̂ = (X₀' * W * X₀)\(X₀' * W * y)
    return β̂
end

function local_poly_reg!(X₀, x̄, W, x, y, x₀, h)
    copyto!(x̄, x .- x₀)
    X₀ = build_X₀!(X₀, x̄, 2)
    K(u) = pdf(Epanechnikov(), u)
    copyto!(W, (w = Diagonal(K.(x̄ ./ h)); w/sum(w)))
    β̂ = (X₀' * W * X₀)\(X₀' * W * y)
    return β̂
end

nearest_neighbors(x̄, N) = sort(x̄, by=abs)[1:N]

end
