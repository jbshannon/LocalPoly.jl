module LocalPoly

using Distributions
using LinearAlgebra
using LoopVectorization
using Parameters
using StaticArrays
using StatsBase

import Polynomials
import Base.show

export LPModel
export lpreg

include("binning.jl")
include("kernels.jl")
include("regressions.jl")

end
