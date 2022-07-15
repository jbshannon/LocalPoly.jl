"""
$(README)

---
## Exports
$(EXPORTS)
"""
module LocalPoly

using Distributions
using DocStringExtensions
using LinearAlgebra
using LoopVectorization
using Parameters
using StaticArrays
using StatsBase

import Polynomials
import Base.show

export LPModel
export lpreg!, lpreg
export plugin_bandwidth

include("binning.jl")
include("kernels.jl")
include("regressions.jl")

end
