"`Dict` containing supported kernel functions"
const KERNELS = Dict(
    :Uniform => u -> abs(u) <= 1 ? 0.5 : 0.0,
    :Triangular => u -> abs(u) <= 1 ? 1-abs(u) : 0.0,
    :Epanechnikov => u -> abs(u) <= 1 ? 3*(1-u^2)/4 : 0.0,
    :Quartic => u -> abs(u) <= 1 ? 15*((1-u^2)^2)/16 : 0.0,
    :Triweight => u -> abs(u) <= 1 ? 35*((1-u^2)^3)/32 : 0.0,
    :Tricube => u -> abs(u) <= 1 ? 70*((1-abs(u)^3)^3)/81 : 0.0,
    :Gaussian => u -> pdf(Normal(), u),
    :Cosine => u -> abs(u) <= 1 ? (π/4)*cos((π/2)*u) : 0.0,
    :Logistic => u -> 1/(exp(u) + 2 + exp(-u)),
    :Parzen => u -> (ū=abs(u); ū <= 1 ? ū <= 0.5 ? 4/3-8ū^2+8ū^3 : 8*(1-ū^3)/3 : 0.0),
    :Sigmoid => u -> (2/π)/(exp(u) + exp(-u)),
    :Silverman => u -> (ū = abs(u)/√2; 0.5*exp(-ū)*sin(ū + π/4)),
)

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

    (0, 1, :Biweight) => 2.036,
    (0, 3, :Biweight) => 3.633,
    (1, 2, :Biweight) => 2.586,
    (2, 3, :Biweight) => 3.208,

    (0, 1, :Triweight) => 2.312,
    (0, 3, :Triweight) => 3.987,
    (1, 2, :Triweight) => 2.869,
    (2, 3, :Triweight) => 3.503,
)

"Fan-Gijbels plugin bandwidth estimator"
function plugin_bandwidth(
    x::AbstractVector, y::AbstractVector, ν::Int, p::Int;
    kernel=:Epanechnikov
)
    m̌ = Polynomials.fit(x, y, p+3)
    m̌⁽ᵖ⁺¹⁾ = Polynomials.derivative(m̌, p+1)
    ε̂ = @turbo @. y - m̌(x)
    σ̃² = var(ε̂)
    ȟ = 𝐶[(ν, p, kernel)] * (σ̃² / sum(abs2∘m̌⁽ᵖ⁺¹⁾, x))^(1/(2p+3))
    return ȟ
end
