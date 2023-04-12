"""
$(README)

---
## Exports
$(EXPORTS)
"""
module LocalPoly

using Distributions
using DocStringExtensions
using IfElse
using LinearAlgebra
using LoopVectorization
using Parameters
using StaticArrays
using StatsBase

import Base.show

export LPModel, LPGridModel
export lpreg!, lpreg
export plugin_bandwidth

include("binning.jl")
include("kernels.jl")
include("regressions.jl")
include("gridmodel.jl")

end
