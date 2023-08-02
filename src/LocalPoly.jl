"""
$(README)

---
## Exports
$(EXPORTS)
"""
module LocalPoly

using Base.Cartesian
using Combinatorics
using Distributions
using DocStringExtensions
using IfElse
using FFTW
using LinearAlgebra
using LoopVectorization
using Parameters
using StaticArrays
using StatsBase

import Base.show

export linear_binning, linear_binning!, gridnodes
export plugin_bandwidth
export lpreg!, lpreg
export confint

include("binning.jl")
include("kernels.jl")
include("polynomials.jl")
include("regression.jl")
include("convolution.jl")
include("fft.jl")
include("MSE.jl")

function lpreg(
    y, x;
    degree=1,
    nbins=guessbins(x),
    kernel=:Epanechnikov,
    h=plugin_bandwidth(linear_binning(x, y; nbins); ν=degree-1, p=degree, kernel),
)
    grid = linear_binning(x, y; nbins)
    convdata = ConvolutionData(grid, h; degree, kernel)
    fftdata = FFTData(convdata)
    convolve!(fftdata)
    β̂ = finalstep(convdata)
    return β̂, gridnodes(grid)
end

function lpreg(
    grid::GridData{T, N, R};
    degree=1,
    kernel=:Epanechnikov,
    h=plugin_bandwidth(grid; ν=degree-1, p=degree, kernel),
) where {T, N, R}
    convdata = ConvolutionData(grid, h; degree, kernel)
    fftdata = FFTData(convdata)
    convolve!(fftdata)
    β̂ = finalstep(convdata)
    return β̂
end

function lpreg(
    y, x, v;
    degree=1,
    nbins=0,
    kernel=:Epanechnikov,
    h=plugin_bandwidth(linear_binning(x, y; nbins); ν=degree-1, p=degree, kernel),
)
    if nbins > 0
        grid = linear_binning(x, y; nbins)
        𝐑 = RegressionData(grid, degree)
    else
        𝐑 = RegressionData(y, x, degree)
    end
    β̂ = lpreg!(𝐑, v; kernel, h)
    return β̂
end

end # module
