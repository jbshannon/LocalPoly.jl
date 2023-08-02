@testset verbose=true "Kernels" begin
    x = -10:0.0001:10
    kernels = [
        :Uniform,
        :Triangular,
        :Epanechnikov,
        :Quartic,
        :Triweight,
        :Tricube,
        :Gaussian,
        :Cosine,
        :Logistic,
        :Sigmoid,
        :Silverman,
    ]
    for kname in kernels
        K(u) = LocalPoly.kernelfunc(Val(kname), u)
        # Test kernel function is mean 0
        @test sum(u -> u*K(u), x) ≈ 0 atol=1e-4

        # Test kernel function integrates to 1
        @test sum(u -> x.step.hi*K(u), x) ≈ 1 atol=1e-2
    end
end
