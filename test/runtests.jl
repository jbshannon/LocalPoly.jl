using LocalPoly
using StatsBase
using Test

@testset verbose=true "LocalPoly.jl" begin

    # @testset verbose=true "Design matrix" begin
    #     # Check that design matrix function works as planned
    #     y = zeros(3)
    #     x = [1, 2, 3]
    #     M = LocalPoly.RegressionData(y, x, 1)
    #     X = LocalPoly._polybasis(x, 0, 2)
    #     @test all(X .== hcat(x .^0, x, x .^2))

    #     LocalPoly._polybasis!(X, x, 1)
    #     x̂ = x .- 1
    #     @test all(X .== hcat(x̂ .^0, x̂, x̂ .^2))
    # end

    include("kernels.jl")
    include("binning.jl")
    include("regressions.jl")
end
