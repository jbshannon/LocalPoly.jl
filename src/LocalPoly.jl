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
using Reexport
using StatsBase

@reexport using StaticArrays

import Polynomials
import Base.show

export LPModel
export lpreg!, lpreg
export plugin_bandwidth

include("binning.jl")
include("kernels.jl")
include("regressions.jl")

end
